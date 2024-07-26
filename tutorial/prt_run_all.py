import pyrotini as prt

def run_all():
    # Parsing 
    args = prt.prt_parser()
    pdb = args.pdb
    wdir = args.wdir
    
    # Preparing files
    prt.prepare_files(pdb, wdir)
    
    # Coarse-graining the protein + solvating
    prt.martinize_go('protein.pdb', wdir, go_map='go.map', go_eps=10.0, go_moltype="protein", go_low=0.3, go_up=0.8, go_res_dist=3)
    prt.solvate(wdir, bt='dodecahedron', d=1.25, radius=0.21, conc=0.0)
    
    # MD
    prt.energy_minimization(wdir, ncpus=0)
    prt.heatup(wdir, ncpus=0)
    prt.equilibration(wdir, ncpus=0)
    prt.md(wdir, nsteps=-2, ncpus=0)
    
    # Analysis
    prt.convert_trajectory(wdir, tb=0, te=1000) # time in ns 
    prt.get_covariance_matrix(wdir)
    prt.get_dfis(wdir)
    prt.plot_dfi(wdir)
  
    
if __name__ == '__main__':
    run_all()