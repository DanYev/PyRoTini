import os
import pyrotini as prt

def main():
    args = prt.prt_design_parser()
    pdb_file = args.pdb
    target = args.target
    run = args.run

    print(f"Starting design for {pdb_file}, targeting residue {target}, run {run}")

    prt.initialize_pyrosetta()
    original_pose, starting_pose, name = prt.clean_and_load_pdb(pdb_file)
    scorefxn = prt.pyrosetta.get_fa_scorefxn()
    print(f"Original score: {scorefxn(original_pose):.2f}")
    
    main_dir, output_dir = prt.setup_dir(target)
    relaxed_pose = prt.fast_relax(starting_pose, scorefxn)
    
    job = prt.setup_job_distributor(f'{name}_design_{run}', 10, scorefxn)
    job.native_pose = original_pose
    pose = prt.pyrosetta.Pose()

    print("Starting design process...")
    while not job.job_complete:
        pose.assign(relaxed_pose)
        prt.design_around(pose, target, scorefxn)
        mutations = prt.print_mutations(original_pose, pose, job.current_name)
        if not mutations:
            continue
        prt.minimize(pose, scorefxn)
        prt.fast_relax(pose, scorefxn)
        job.output_decoy(pose)
    os.chdir(main_dir)
    print(f"PyRosetta design completed. Output files are in {output_dir}")

if __name__ == '__main__':
    main()

# python prt_design.py -f 1btl.pdb -t 44A -r 1