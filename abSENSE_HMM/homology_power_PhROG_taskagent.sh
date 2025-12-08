#!/bin/sh

for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
  for LOO_HMMSEARCH_FILENAME in $hmm_LOO_results_prefix/results/hmmsearch_out/${OGID}*_LOO_hmmsearch.out
  do
    PHROGID=${LOO_HMMSEARCH_FILENAME%_LOO_hmmsearch.out}
    PHROGID=${PHROGID##*/}
    Rscript run_absense.R $PHROGID $hmm_LOO_results_prefix $alerax_dir $posterior_clades_dir $outdir $BOOL_DELINEATE_PARALOGS
  done
done;
