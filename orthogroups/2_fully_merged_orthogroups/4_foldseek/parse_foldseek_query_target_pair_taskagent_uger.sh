#!/bin/sh

for QUERY_OG in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`
do
    for TARGET_OG in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 2`
    do
        INPUTFILE="$DATADIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10.tsv"
        OUTFILE="$OUTDIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10_mean.bitscore.tsv"

        if [[ -s "$INPUTFILE" && ! -f "$OUTFILE" ]]
        then
            python parse_foldseek_query_target_pair.py $INPUTFILE $OUTFILE $PDBDIR
        fi
    done
done
