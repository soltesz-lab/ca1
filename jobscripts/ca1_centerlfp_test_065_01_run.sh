#!/bin/bash
### set the number of nodes and the number of PEs per node
#PBS -q normal
#PBS -l nodes=512:ppn=16:xe
### set the wallclock time
#PBS -l walltime=2:00:00
### set the job name
#PBS -N ca1_centerlfp_test_065_01
### set the job stdout and stderr
#PBS -e ./results/$PBS_JOBID.err
#PBS -o ./results/$PBS_JOBID.out
### set email notification
### Set umask so users in my group can read job stdout and stderr files
#PBS -W umask=0027
#PBS -A bayj

module swap PrgEnv-cray PrgEnv-gnu
module load bwpy 
module load bwpy-mpi

set -x

## go to directory from where job script was launched
cd $PBS_O_WORKDIR

export SCRATCH=/projects/sciteam/bayj
export NEURONROOT=$SCRATCH/nrn
export PYTHONPATH=$NEURONROOT/lib/python:$PYTHONPATH
export PATH=$NEURONROOT/x86_64/bin:$PATH

aprun -n 8192 -b -- bwpy-environ -- \
nrniv -mpi jobscripts/ca1_centerlfp_test_065_01_run.hoc -c "quit()"
cp ./jobscripts/ca1_centerlfp_test_065_01* ./results/ca1_centerlfp_test_065_01

