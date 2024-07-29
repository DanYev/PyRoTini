import argparse
import importlib.resources
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import pyrosetta
import pandas as pd
import shutil
import subprocess as sp
from pathlib import Path
from pyrotini.get_go import get_go
from . import PRT_DATA, PRT_DICT


def prt_parser():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Runs CGMD simulation for pyRosetta")
    parser.add_argument("-f", "--pdb", required=True, help="PDB file")
    parser.add_argument("-d", "--wdir", required=True, help="Relative path to the working directory")
    parser.add_argument("-n", "--ncpus", required=False, help="Number of requested CPUS for a batch job")
    return parser.parse_args()

def prt_design_parser():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Runs PyRosetta design")
    parser.add_argument("-f", "--pdb", required=True, help="PDB file")
    parser.add_argument("-t", "--target", type=str, required=True, help="Target residue to mutate around as integer")
    parser.add_argument("-r", "--run", type=int, required=True, help="Run number for the design")
    return parser.parse_args()

def initialize_pyrosetta():
    """Initialize PyRosetta with custom flags."""
    pyrosetta.init('''
        -relax:default_repeats 5
        -relax:constrain_relax_to_start_coords
        -relax:coord_constrain_sidechains
        -relax:ramp_constraints false
        -score::hbond_params correct_params
        -no_his_his_pairE
        -extrachi_cutoff 1
        -multi_cool_annealer 10
        -ex1 -ex2
        -use_input_sc
        -flip_HNQ
        -ignore_unrecognized_res
        -relax:coord_cst_stdev 0.5
    ''')

def clean_and_load_pdb(filename):
    """Clean and load the PDB file."""
    pyrosetta.toolbox.cleaning.cleanATOM(filename)
    name = os.path.splitext(filename)[0]
    starting_pose = pyrosetta.pose_from_pdb(f"{name}.clean.pdb")
    original_pose = starting_pose.clone()
    
    return original_pose, starting_pose, name

def setup_dir(target):
    """Set up the design environment and return necessary objects."""
    main_dir = os.getcwd()
    output_dir = os.path.join("Outputs", f"mutating_around_res{target}")
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)
    
    return main_dir, output_dir 

def fast_relax(pose, scorefxn):
    fr = pyrosetta.rosetta.protocols.relax.FastRelax()
    fr.set_scorefxn(scorefxn)
    fr.apply(pose)
    print(f"Relaxed score: {scorefxn(pose):.2f}")
    relaxed_pose = pose.clone()

    return relaxed_pose

def setup_job_distributor(name, n_decoys, scorefxn):
    """Set up the job distributor."""
    return pyrosetta.toolbox.py_jobdistributor.PyJobDistributor(name, n_decoys, scorefxn)

def design_around(pose, target, scorefxn):
    """Set up and apply the design protocol."""
    target_residue = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(str(target))
    design_radius = random.randint(8, 12)
    design_shell = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(target_residue, design_radius, False)
    repack_shell = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(target_residue, design_radius + 4, True)
    imp_residues = "2,19,20,37,42,45,48,51,82,97,105,107,120,141,142,149,158,194,201,209,218,226,230,235"
    imp_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(imp_residues)

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), imp_residue_selector, False))
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), design_shell, True))
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), repack_shell, True))

    prm = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    prm.task_factory(tf)
    prm.score_function(scorefxn)
    prm.apply(pose)
    
    print(f"Design: {scorefxn(pose):.2f}")
    print_radius(design_radius)

def minimize(pose, scorefxn):
    """Set up the minimize movers and minimize pose for the design process."""
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)
    mm.set_jump(True)

    min_bb = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    min_bb.score_function(scorefxn)
    min_bb.movemap(mm)
    min_bb.min_type('dfpmin_armijo_nonmonotone')
    min_bb.tolerance(0.001)
    min_bb.max_iter(1000)

    multi_min = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()
    multi_min.set_mover(min_bb)
    multi_min.set_scorefxn(scorefxn)
    multi_min.set_max_accepted_trials(10)
    multi_min.set_sampletype("low")
    multi_min.set_temperature(0.6)
    multi_min.set_recover_low(True)
    multi_min.set_preapply(False)

    multi_min.apply(pose)
    print(f"After minimization score: {scorefxn(pose):.2f}")

def print_mutations(original_pose, designed_pose, job_name):
    """Print and log the mutations between the original and designed poses."""
    original_seq = original_pose.sequence()
    designed_seq = designed_pose.sequence()
    mutations = []
    for i in range(0, original_pose.total_residue()):
        if original_seq[i] != designed_seq[i]:
            resid = original_pose.pdb_info().pose2pdb(i+1)
            mutations.append(f"{original_seq[i]}{resid.split()[0]}{designed_seq[i]}")
    mutations_str = f"{job_name} - Mutations: " + ", ".join(mutations)
    print(mutations_str)
    with open("mutations.txt", "a") as mutation_file:
        if mutations:
            mutation_file.write(mutations_str + "\n")
        else:
            mutation_file.write("\n")
    return mutations

def print_radius(radius):
    """Print and log the radius of design shell and repack shell."""
    radius_str = f"Design shell: {radius}A      Repack shell: {radius+4}A"
    print(radius_str)
    with open("radius.txt", "a") as radius_file:
        radius_file.write(radius_str + "\n")
                  
def make_topology_file(wdir, protein='protein'):
    r"""
    -protein        Name of the protein (just for the reference, doesn't affect anything)
    
    """
    bdir = os.getcwd()
    os.chdir(wdir)
    with open('system.top', 'w') as out_file:
        out_file.write(f'#define GO_VIRT\n')
        out_file.write(f'#include "martini.itp"\n')
        out_file.write(f'#include "go_atomtypes.itp"\n')
        out_file.write(f'#include "go_nbparams.itp"\n')
        out_file.write(f'#include "protein.itp"\n')
        out_file.write(f'#include "solvents.itp"\n')
        out_file.write(f'#include "ions.itp"\n')
        out_file.write(f'\n[ system ]\n')
        out_file.write(f'Martini protein in water\n\n') 
        out_file.write(f'\n[ molecules ]\n')
        out_file.write(f'{protein}  1\n')
    os.chdir(bdir)
        
        
def link_itps(wdir):
    bdir = os.getcwd()
    os.chdir(wdir)
    for name, path in PRT_DICT.items():
        if name.endswith('.itp'):
            command = f'ln -sf {path} {name}' # > /dev/null 2>&
            sp.run(command.split())
    os.chdir(bdir)
    
    
def fix_go_map(wdir, in_map, out_map='go.map'):
    bdir = os.getcwd()
    os.chdir(wdir)
    with open (in_map, 'r') as in_file:
         with open (out_map, 'w') as out_file:
            for line in in_file:
                if line.startswith('R '):
                    new_line = ' '.join(line.split()[:-1])
                    out_file.write(f'{new_line}\n')
    os.chdir(bdir)


def prepare_files(pdb, wdir='test', protein='protein'):
    r"""
    -wdir           Relative path to the working directory
    """
    os.makedirs(wdir, exist_ok=True)
    copy_to = os.path.join(wdir, 'protein.pdb')
    shutil.copy(pdb, copy_to)
    link_itps(wdir)
    make_topology_file(wdir, protein=protein)
    print("Getting Go-map...")
    get_go(wdir, protein)
    fix_go_map(wdir, in_map='protein.map')
    print('All the files are ready!')
    
    
def martinize_go(pdb, wdir, go_map='go.map', go_eps=10.0, go_moltype="protein", go_low=0.3, go_up=0.8, go_res_dist=3):
    r"""
    Virtual site based GoMartini:
    -go_map         Contact map to be used for the Martini Go model.Currently, only one format is supported. (default: None)
    -go_moltype     Set the name of the molecule when using Virtual Sites GoMartini. (default: protein)
    -go_eps         The strength of the Go model structural bias in kJ/mol. (default: 9.414)                        
    -go_low         Minimum distance (nm) below which contacts are removed. (default: 0.3)
    -go_up          Maximum distance (nm) above which contacts are removed. (default: 1.1)
    -go_res_dist    Minimum graph distance (similar sequence distance) below whichcontacts are removed. (default: 3)
    """
    bdir = os.getcwd()
    os.chdir(wdir)
    shutil.copy(pdb, 'protein_aa.pdb')
    command = f'martinize2 -f {pdb} -go {go_map} -go-moltype {go_moltype} -go-eps {go_eps} \
        -go-low {go_low} -go-up {go_up} -go-res-dist {go_res_dist} \
        -o protein.top -x protein.pdb -p backbone -dssp -ff martini3001 \
        -sep -scfix -cys 0.3 -resid input -maxwarn 1000'
    sp.run(command.split())
    os.chdir(bdir)
    
    
def solvate(wdir, bt='dodecahedron', d=1.25, radius=0.21, conc=0.0):
    r"""
    -bt             Box type for -box and -d: triclinic, cubic, dodecahedron, octahedron (default: triclinic)
    -d              Distance between the solute and the box (default: 1.25 nm)
    -radius         VWD radius (default: 0.21 nm)
    -conc           Ionic concentration (micromol/l) (default: 0.0 nm)
    """
    bdir = os.getcwd()
    os.chdir(wdir)
    command = f'gmx_mpi editconf -f protein.pdb -c -bt {bt} -d {d} -o system.gro'
    sp.run(command.split())
    command = f'gmx_mpi solvate -cp system.gro -cs {PRT_DICT['water.gro']} -p system.top -radius {radius} -o system.gro'
    sp.run(command.split())
    command = f'gmx_mpi grompp -f {PRT_DICT['ions.mdp']} -c system.gro -p system.top -o ions.tpr -maxwarn 1000'
    sp.run(command.split())
    command = f'gmx_mpi genion -s ions.tpr -p system.top -conc {conc} -neutral -pname NA -nname CL -o system.gro'
    sp.run(command.split(), input='W\n', text=True)
    os.chdir(bdir)
    
    
def energy_minimization(wdir, ncpus=0):
    """
    Perform energy minimization using GROMACS.

    Parameters:
    wdir (str): The working directory where the energy minimization will be performed.
    ncpus (int, optional): Number of CPU threads to use for the minimization. Defaults to 0, 
                           which lets GROMACS decide the number of threads.

    Raises:
    FileNotFoundError: If the necessary input files are not found in the specified directories.
    RuntimeError: If the GROMACS commands fail to execute.
    """
    bdir = os.getcwd()
    os.chdir(wdir)
    os.makedirs('mdrun', exist_ok=True)
    os.chdir('mdrun')
    command = f'gmx_mpi grompp -f {PRT_DICT['em.mdp']} -c ../system.gro -r ../system.gro -p ../system.top -o em.tpr'
    sp.run(command.split())
    options = f'-ntomp {ncpus} -pin on -pinstride 1'
    command = f'gmx_mpi mdrun {options} -deffnm em'
    sp.run(command.split())
    os.chdir(bdir)
    
    
def heatup(wdir, ncpus=0):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.chdir('mdrun')
    command = f'gmx_mpi grompp -f {PRT_DICT['hu.mdp']} -c em.gro -r em.gro -p ../system.top -o hu.tpr'
    sp.run(command.split())
    options = f'-ntomp {ncpus} -pin on -pinstride 1'
    command = f'gmx_mpi mdrun {options} -deffnm hu'
    sp.run(command.split())
    os.chdir(bdir)


def equilibration(wdir, ncpus=0):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.chdir('mdrun')
    command = f'gmx_mpi grompp -f {PRT_DICT['eq.mdp']} -c hu.gro -r hu.gro -p ../system.top -o eq.tpr'
    sp.run(command.split())
    options = f'-ntomp {ncpus} -pin on -pinstride 1'
    command = f'gmx_mpi mdrun {options} -deffnm eq'
    sp.run(command.split())
    os.chdir(bdir)
    

def md(wdir, nsteps=-2, ncpus=0):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.chdir('mdrun')
    command = f'gmx_mpi grompp -f {PRT_DICT['md.mdp']} -c eq.gro -r eq.gro -p ../system.top -o md.tpr'
    sp.run(command.split())
    options = f'-ntomp {ncpus} -pin on -pinstride 1'
    command = f'gmx_mpi mdrun {options} -deffnm md -nsteps {nsteps}'
    sp.run(command.split())
    os.chdir(bdir)
    

def convert_trajectory(wdir, tb=0, te=1000):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.makedirs('analysis', exist_ok=True)
    analysis_dir = os.path.abspath('analysis')
    os.chdir('mdrun')
    shutil.copy('md.tpr', f'{analysis_dir}/md.tpr')
    command = f'gmx_mpi trjconv -s md.tpr -f md.xtc -o {analysis_dir}/mdc.xtc -fit rot+trans -b {tb} -e {te} -tu ns'
    sp.run(command.split(), input='1\n1\n', text=True)
    os.chdir(bdir)
    

def get_covariance_matrix(wdir, tb=0, te=250, tw=5, td=5):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.chdir('analysis')
    for t in range(tb, te-tw-td, td):
        b = t
        e = t + tw
        command = f'gmx_mpi covar -s md.tpr -f mdc.xtc -ascii covar_{t}.dat -b {b} -e {e} -tu ns -last 50'
        sp.run(command.split(), input='1\n3\n', text=True)
    files_to_delete = [f for f in os.listdir() if f.startswith('#')]
    for f in files_to_delete:
        os.remove(f)
    os.chdir(bdir)


def parse_covarince_matrix(file):
    print(f"Reading covariance matrix from {file}")
    df = pd.read_csv(file, sep='\\s+', header=None)
    covarince_matrix = np.asarray(df, dtype=np.float64)
    resn = int(np.sqrt(len(covarince_matrix) / 3)) # number of residues
    covarince_matrix = np.reshape(covarince_matrix, (3*resn, 3*resn))
    return covarince_matrix, resn
    
    
def calculate_dfi(covarince_matrix, resn):
    print("Calculating DFI")
    cov = covarince_matrix
    directions = np.array(([1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,1,1]), dtype=np.float64)
    norm = np.sqrt(np.sum(directions, axis=1))
    directions = (directions.T / norm).T
    dfi = np.zeros(resn)
    for f in directions:
        f = np.tile(f, resn)
        m = covarince_matrix * f
        m = np.reshape(m, (resn, resn, 3, 3))
        m = np.sum(m, axis=-1)
        pert_mat = np.sqrt(np.sum(m * m, axis=-1))
        dfi += np.sum(pert_mat, axis=-1)
    dfi /= np.sum(dfi)
    return dfi
 
 
def percentile(x):
    """
    Calculate the percentile ranking of each element in a 1-dimensional array.

    Parameters:
    x (np.ndarray): Input array.

    Returns:
    np.ndarray: Array of percentile rankings.
    """
    sorted_x = np.argsort(x)
    px = np.zeros(len(x))
    for n in range(len(x)):
        px[n] = np.where(sorted_x == n)[0][0] / len(x)
    return px 
            
            
def get_dfi(file, ):
    covarince_matrix, resn = parse_covarince_matrix(file)
    dfi = calculate_dfi(covarince_matrix, resn)
    pdfi = percentile(dfi)
    print("Saving DFI")
    data = pd.DataFrame()
    resnums = list(range(26, resn + 26))
    data["resn"] = resnums
    data["dfi"] = dfi
    data["pdfi"] = pdfi
    data.to_csv(file.replace('covar', 'dfi'), index=False, header=None, sep=' ')
    return resnums, dfi


def get_dfis(wdir):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.chdir('analysis')
    covar_files = [f for f in os.listdir() if f.startswith('covar') and f.endswith('.dat')]
    dfis = []
    for file in covar_files:
        resnums, dfi = get_dfi(file)
        dfis.append(dfi)
    dfi_av = np.average(np.array(dfis).T, axis=-1)
    data = pd.DataFrame()
    data["resn"] = resnums
    data["dfi"] = dfi_av
    data["pdfi"] = percentile(dfi_av)
    data.to_csv('dfi_av.dat', index=False, header=None, sep=' ')
    os.chdir(bdir)
    
    
def parse_dfi(file):
    df = pd.read_csv(file, sep=' ', header=None)
    res = df[0]
    dfi = df[1]
    pdfi = df[2]
    return res, dfi, pdfi
    
    
def plot_dfi(wdir):
    bdir = os.getcwd()
    os.chdir(wdir)
    os.chdir('analysis')
    dfi_files = [f for f in os.listdir() if f.startswith('dfi_av') and f.endswith('.dat')]
    fig = plt.figure(figsize=(12,4))
    for file in dfi_files:
        res, dfi, pdfi = parse_dfi(file)
        plt.plot(res, pdfi)
    plt.autoscale(tight=True)
    plt.grid()
    plt.tight_layout()
    fig.savefig(f'dfi.png')
    plt.close()
    os.chdir(bdir)    
    

if __name__  == "__main__":
    pass



