#!/bin/bash

for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	FAAFILE="$FAADIR/$OGID.faa"
	MSAFILE="${FAAFILE}.msa"
	
	# Align
	conda activate mafft
	mafft --auto --quiet --thread 1 --anysymbol $FAAFILE > $MSAFILE

	# Trim
	conda activate clipkit_2.3.0
	clipkit -q -m gappy "${MSAFILE}" -o "${FAAFILE}_clipkit.gappy.msa"
	python utils/remove_gap_only_seqs.py "${FAAFILE}_clipkit.gappy.msa"
done;
