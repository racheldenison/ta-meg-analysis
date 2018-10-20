#!/bin/bash

#SBATCH --job-name=meg_stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5GB
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rnd224@nyu.edu
#SBATCH --output=slurm_%A_%a.out
 
module purge
module load r/intel/3.4.2
 
cd ~/code/TA_MEG
R --no-save -q -f meg_timeseriesPermutationTests_hpc.R