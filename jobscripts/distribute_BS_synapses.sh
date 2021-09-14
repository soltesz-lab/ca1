#!/bin/bash
mpirun.mpich -np 1 python3 ./scripts/distribute_synapse_locs.py \
             --template-path templates \
              --config=Full_Scale.yaml \
              --populations BS \
              --forest-path=./datasets/BS_forest.h5 \
              --output-path=./datasets/BS_forest_syns_20210912.h5 \
              --distribution=poisson \
              --io-size=1 -v
