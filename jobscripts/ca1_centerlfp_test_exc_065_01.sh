#!/bin/bash
#SBATCH --job-name="ca1_centerlfp_test_exc_065_01"
#SBATCH --output="./jobscripts/ca1_centerlfp_test_exc_065_01.%N.o"
#SBATCH --partition=normal
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=56
#SBATCH -p normal
#SBATCH -t 2:00:00


set -x

mkdir -p ./results/ca1_centerlfp_test_exc_065_01/

ibrun x86_64/special -mpi ./jobscripts/ca1_centerlfp_test_exc_065_01_run.hoc

cp ./jobscripts/ca1_centerlfp_test_exc_065_01* ./results/ca1_centerlfp_test_exc_065_01/
