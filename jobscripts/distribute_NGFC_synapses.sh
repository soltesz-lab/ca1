#!/bin/bash
mpirun.mpich -np 1 python3 ./scripts/distribute_synapse_locs.py \
             --template-path templates \
              --config=Full_Scale.yaml \
              --populations NGFC \
              --forest-path=./datasets/NGFC_forest.h5 \
              --output-path=./datasets/NGFC_forest_syns_20211026.h5 \
              --distribution=poisson \
              --io-size=1 -v
