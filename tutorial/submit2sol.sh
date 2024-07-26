#!/bin/bash
#SBATCH --time=0-01:00:00                                                       # upper bound time limit for job to finish d-hh:mm:ss
#SBATCH --partition=general
#SBATCH --qos=public                                                            # public grp_sozkan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --gres=gpu:1	                                                        # number of GPUs
#SBATCH -o slurm_output/output.%A.%a.out
#SBATCH -e slurm_output/error.%A.%a.err

python prt_run_all.py -f '1btl.pdb' -d 'test' 
