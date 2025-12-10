#!/bin/sh

for TAXID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	fasta-36.3.8i/bin/ggsearch36 -E 1 -m 8 -T 2 $QUERY_FASTA $TARGET_FASTA > $OUTPUT_FILE
done;
