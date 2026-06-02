#!/bin/bash

conda activate graphpart
source graph-part/.venv/bin/activate

FAAFILE="zenodo_directory/data/deeploc/swissprot_PCP_All_2026.04.05.fasta"
OUTDIR="output_directory"

# Precomputed, ignore priority
graphpart precomputed --fasta-file $FAAFILE --threshold 0.3 --out-file $OUTDIR/graphpart_assignments.csv --labels-name label --partitions 5 --edge-file data/deeploc/ggsearch/ggsearch_swissprot_mitoepi_combined_fasta_vs_self_expect1_pident0.3.csv --metric-column 2

