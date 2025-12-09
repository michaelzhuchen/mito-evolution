#!/bin/bash

for OG_ID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	# Run MMseqs2
	mmseqs easy-cluster ${OG_ID}.faa ${OUTPUT_DIRECTORY} ${TEMPORARY_DIRECTORY} --threads 1 --split-memory-limit 2G --min-seq-id 0.5 --cov-mode 2 -c 0.8 -v 1
done;
