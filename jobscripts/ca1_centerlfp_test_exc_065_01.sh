#!/bin/bash
#SBATCH --job-name="ca1_centerlfp_test_exc_065_01"
#SBATCH --output="./jobscripts/ca1_centerlfp_test_exc_065_01.%N.o"
#SBATCH --partition=compute
#SBATCH --nodes 4
#SBATCH --ntasks-per-node=12
#SBATCH --export=ALL
#SBATCH -t 0:30:00

module load python
module load mpi4py

set -x

ibrun -v x86_64/special -mpi ./jobscripts/ca1_centerlfp_test_exc_065_01_run.hoc
cp ./jobscripts/ca1_centerlfp_test_exc_065_01* ./results/ca1_centerlfp_test_exc_065_01/
