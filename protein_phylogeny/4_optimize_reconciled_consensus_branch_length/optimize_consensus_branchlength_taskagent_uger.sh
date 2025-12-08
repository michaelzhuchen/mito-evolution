#!/bin/sh

reconciled_consensus_tree_directory="data/reconciled_trees"

for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	MSAFILE="${OGID}.faa_clipkit.gappy.msa"

	CONSTRAINT_TREE_FILE="${reconciled_consensus_tree_directory}/$OGID/reconciliations/summaries/family_1_consensus_50.newick"

	OUTPUT_PREFIX="$OUTDIR/${OGID}_consensus50_blopt"
	IQTREE_OUTPUT="${OUTPUT_PREFIX}.iqtree"

	ORIGINAL_IQTREE_FILE="$ORIGINALDIR/${MSAFILE}.iqtree"

	# Get best fit substitution model
	SUBSTITUTION_MODEL=$(sed -n 's/^Best-fit model according to BIC: //p' $ORIGINAL_IQTREE_FILE)

	# Prune MSA to remove sequences that were pruned during infer origin rooting
	Rscript prune_msa.R $ORIGINALDIR $ORIGINROOTDIR $OUTDIR $OGID
	PRUNEDMSAFILE="$OUTDIR/${OGID}.faa_clipkit.gappy_pruned.msa"

	## Constraint tree (resolve polytomies)
	iqtree2 -s $PRUNEDMSAFILE -g $CONSTRAINT_TREE_FILE --prefix $OUTPUT_PREFIX -m $SUBSTITUTION_MODEL -nt $N_CORES -seed 42 -quiet
	
done;
