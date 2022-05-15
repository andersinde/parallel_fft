#!/bin/bash -l

# time allocation
#SBATCH -A edu22.sf2568
#SBATCH -J anders_fft
#SBATCH -t 00:01:00


#SBATCH --ntasks=64
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32


# Number of tasks per core (prevent hyperthreading)
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
# Use Dardel's main partition
#SBATCH -p main

#srun ./program 11
#srun ./program 12
#srun ./program 13
#srun ./program 14
#srun ./program 15
#srun ./program 16
#srun ./program 17
#srun ./program 18
srun ./program 19
srun ./program 20 # segmentation fault :(
