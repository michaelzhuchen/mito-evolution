#!/bin/bash

DIR="ACANB_virus_alignments_and_trees"

OUTFILE="${DIR}_highest_LL_trees.txt"
rm -rf $OUTFILE

OUTDIR="${DIR}_highest_LL_trees"
mkdir -p $OUTDIR

# Find all completed treefile prefixes using .iqtree files
treefile_list="${DIR}_completed_treefile_prefixes.txt"
find $DIR -type f -name "*.iqtree" -exec sh -c 'echo "$(basename "${1%.iqtree}")"' _ {} \; > $treefile_list

while read line
do
	MSAFILE="$DIR/$line"
	IQTREE_LOGFILE="$MSAFILE.iqtree"
    # If IQTREE job finished successfully
	if grep -q "Total wall-clock time used" "$IQTREE_LOGFILE"; then
        
        # If the consensus tree had a higher LL than the ML tree
		if grep -q "Consensus tree has higher likelihood than ML tree found" "$IQTREE_LOGFILE"; then
			echo "$line.contree" >> $OUTFILE
			cp "$MSAFILE.contree" $OUTDIR
		else
			echo "$line.treefile" >> $OUTFILE
			cp "$MSAFILE.treefile" $OUTDIR
		fi

	fi
done < $treefile_list
