#!/bin/bash

N_THREADS="4"

iqtree2 -s ${MSA_FILE} -m MFP -cmin 4 -cmax 12 -nt ${N_THREADS} -bb 1000 -alrt 1000 -seed 42 -wbtl -bnni -quiet
