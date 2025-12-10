#!/bin/bash

for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	FAAFILE="$FAADIR/$OGID.faa"
	MSAFILE="${FAAFILE}.msa"
	if [[ ! -s $MSAFILE ]]
	then
		SEQ_COUNT=$(grep -c '^>' "$FAAFILE")
		if [ "$SEQ_COUNT" -eq 0 ]; then
			echo "No sequences in $FAAFILE"
		elif [ "$SEQ_COUNT" -le 500 ]; then
			## PPP: most accurate. for <=500 hundred taxa
    		muscle5.1.linux_intel64 -align $FAAFILE -output $MSAFILE -threads $N_CORES -quiet
    	else
			## SUPER5: less accurate, divide and conquer heuristic. for >500 hundred taxa
    		muscle5.1.linux_intel64 -super5 $FAAFILE -output $MSAFILE -threads $N_CORES -quiet
    	fi
    fi
	
	# Trim with clipkit
	TRIMOUTFILE="${FAAFILE}_clipkit.gappy.msa"
	if [[ ! -s $TRIMOUTFILE ]]
	then
		conda activate clipkit_2.3.0
		clipkit -q -m gappy "${MSAFILE}" -o $TRIMOUTFILE
		python utils/remove_gap_only_seqs.py $TRIMOUTFILE
	fi
done;
