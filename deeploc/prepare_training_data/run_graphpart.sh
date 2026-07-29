#!/bin/bash

FAAFILE="mito-evolution/data/deeploc/swissprot_PCP_All_2026.04.05.fasta"
OUTDIR="mito-evolution/data/deeploc/graphpart/swissprot_PCP_All_2026.04.05.fasta_pident0.3_expect1_nopriority"

graphpart precomputed --fasta-file $FAAFILE --threshold 0.3 --out-file $OUTDIR/graphpart_assignments.csv --labels-name label --partitions 5 --edge-file data/deeploc/ggsearch/ggsearch_swissprot_mitoepi_combined_fasta_vs_self_expect1_pident0.3.csv --metric-column 2
