#!/bin/bash
#SBATCH --time=0-04:00:00                                                       # upper bound time limit for job to finish d-hh:mm:ss
#SBATCH --partition=htc
#SBATCH --qos=public                                                            # public grp_sozkan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=gpu:1	                                                        # number of GPUs
#SBATCH -o slurm_output/output.%A.out
#SBATCH -e slurm_output/error.%A.err

python prt_run_all.py -f '1btl.pdb' -d 'test' 
