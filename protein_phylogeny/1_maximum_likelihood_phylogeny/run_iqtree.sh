#!/bin/bash

OGID="MOG0001047"

MSA_FILE="alignments_and_initial_trees/$OGID.faa.msa"
N_THREADS="2"

iqtree2 -s ${MSA_FILE} -m MFP -cmin 4 -cmax 12 -nt ${N_THREADS} -bb 1000 -alrt 1000 -seed 42 -wbtl -bnni -quiet
