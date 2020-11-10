#!/bin/bash
#SBATCH -J ca1_ripplestim_1x_01
#SBATCH -o ./jobscripts/ca1_ripplestim_1x_01.%j.o
#SBATCH -n 3008
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
set -x
ibrun tacc_affinity x86_64/special -mpi ./jobscripts/ca1_ripplestim_1x_01_run.hoc
cp ./jobscripts/ca1_ripplestim_1x_01* ./results/ca1_ripplestim_1x_01/
