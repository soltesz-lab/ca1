#!/bin/bash
mpirun.mpich -np 1 python3 ./scripts/distribute_synapse_locs.py \
             --template-path templates \
              --config=Full_Scale.yaml \
              --populations IVY \
              --forest-path=./datasets/IVY_forest.h5 \
              --output-path=./datasets/IVY_forest_syns_20210925.h5 \
              --distribution=poisson \
              --io-size=1 -v
