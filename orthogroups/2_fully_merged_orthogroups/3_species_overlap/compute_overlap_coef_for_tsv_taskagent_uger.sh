#!/bin/bash

SUFFIX=".faa_clipkit.gappy.msa_global.hhr_prob30"
BOOL_STRUCTURE_MERGE="TRUE"

for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	# Compute the overlap coefficient
	Rscript orthogroups/helpers/compute_overlap_coef_for_tsv.R "$DATADIR" "$OUTDIR" "$OGID" "$SUFFIX" "$BOOL_STRUCTURE_MERGE"
done
