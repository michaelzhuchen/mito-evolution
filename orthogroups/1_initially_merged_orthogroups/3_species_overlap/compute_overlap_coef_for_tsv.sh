#!/bin/bash

DATADIR="output_directory"
OUTDIR="output_directory"
SUFFIX=".faa_clipkit.gappy.msa_global.hhr_expect1e-3"
BOOL_STRUCTURE_MERGE="FALSE"

OGID="OG0002514"

# Compute the overlap coefficient
Rscript orthogroups/helpers/compute_overlap_coef_for_tsv.R "$DATADIR" "$OUTDIR" "$OGID" "$SUFFIX" "$BOOL_STRUCTURE_MERGE"