#!/bin/bash

# Prepare DIAMOND databases and generate list of DIAMOND commands
conda activate OrthoFinder_v255
orthofinder -f $FASTA_DIRECTORY -op

