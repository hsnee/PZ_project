#!/bin/bash -l
#SBATCH --image=docker:lsstdesc/stack-jupyer:weekly-latest
#SBATCH -p regular   #Submit to the regular 'partition'
#SBATCH -N 1     #Use 1 node
#SBATCH -t 00:30:00  #Set up to 30 min time limit
#SBATCH -C knl   #Use KNL nodes
#SBATCH --output=/global/homes/h/husni/mcmc.out
#SBATCH --job-name=mcmc

export OMP_NUM_THREADS=1
srun -n 68 shifter /global/homes/h/husni/PZ_Project/shifter_mcmc.sh
