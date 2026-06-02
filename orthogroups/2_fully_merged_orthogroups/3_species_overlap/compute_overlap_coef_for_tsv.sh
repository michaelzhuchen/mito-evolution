#!/bin/bash

DATADIR="output_directory"
OUTDIR="output_directory"
SUFFIX=".faa_clipkit.gappy.msa_global.hhr_prob30"
BOOL_STRUCTURE_MERGE="TRUE"

OGID="MOG0030340"

# Compute the overlap coefficient
Rscript orthogroups/helpers/compute_overlap_coef_for_tsv.R "$DATADIR" "$OUTDIR" "$OGID" "$SUFFIX" "$BOOL_STRUCTURE_MERGE"