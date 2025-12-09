#!/bin/bash

SUFFIX=".faa_clipkit.gappy.msa_global.hhr_expect1e-3"
BOOL_STRUCTURE_MERGE="FALSE"

for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	# Compute the overlap coefficient
	Rscript orthogroups/helpers/compute_overlap_coef_for_tsv.R "$DATADIR" "$OUTDIR" "$OGID" "$SUFFIX" "$BOOL_STRUCTURE_MERGE"
done
