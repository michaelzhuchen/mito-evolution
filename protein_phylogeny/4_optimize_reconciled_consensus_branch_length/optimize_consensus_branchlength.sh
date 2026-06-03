#!/bin/sh

CONSENSUSDIR="reconciled_trees_species.tree.1"
ORIGINALDIR="alignments_and_initial_trees"
ALERAXDIR="reconciled_trees_species.tree.1"
OUTDIR="reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.1"

OGID="MOG0001047"


MSAFILE="${OGID}.faa_clipkit.gappy.msa"

CONSTRAINT_TREE_FILE="${CONSENSUSDIR}/$OGID/reconciliations/summaries/family_1_consensus_50.newick"

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

## Get clade supports
Rscript get_consensus_clade_supports.R $OUTDIR $ALERAXDIR $OUTDIR $OGID
