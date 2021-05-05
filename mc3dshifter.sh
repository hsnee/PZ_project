#!/bin/bash -l
#SBATCH --image=docker:lsstdesc/stack-jupyter:weekly-latest
#SBATCH -p regular   #Submit to the regular 'partition'
#SBATCH -N 1     #Use 1 node
#SBATCH -t 20:00:00  #Set up to 30 min time limit
#SBATCH -C knl   #Use KNL nodes
#SBATCH --output=/global/u2/h/husni/PZ_Project/shifter-%j.out
#SBATCH --job-name=mcmc3d
#SBATCH --signal=B:USR1@120
export OMP_NUM_THREADS=1
srun -n 64 shifter /global/u2/h/husni/PZ_Project/runshift.sh
