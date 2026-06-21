#!/bin/bash
#SBATCH -N 1            # number of nodes
#SBATCH -c 16           # number of cores 
#SBATCH -p general      # partition
#SBATCH -t 0-06:00:00   # time in d-hh:mm:ss 
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --export=NONE   # Purge the job-submitting shell environment

prot=$1

/packages/envs/pyrosetta-2023/bin/python prt_design.py -t $2 -f ${prot}.pdb -r $3
