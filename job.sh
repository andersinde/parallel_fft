#!/bin/bash -l

# time allocation
#SBATCH -A edu22.sf2568

# #SBATCH --reservation=lab-31-1

# job name
#SBATCH -J anders_fft

# 1 minute wall-clock time will be given to this job
#SBATCH -t 00:01:00

# Number of tasks 
#SBATCH --ntasks=32

#SBATCH --cpus-per-task=1

# Use Dardel's shared partition
#SBATCH -p main

srun ./program 15
srun ./program 16
srun ./program 17
srun ./program 18
srun ./program 19
