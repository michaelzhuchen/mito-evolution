#!/bin/sh

DATASETNAME="species_tree_1"
ALERAX_DIR="reconciled_trees_species.tree.1"
POSTERIOR_CLADES_DIR="reconciled_trees_posterior_clades_species.tree.1"
HMM_LOO_RESULTS_PREFIX="data/abSENSE_HMM/$DATASETNAME"
SPECIES_TREE_FILENAME="data/species_phylogeny/processed_species_tree/${DATASETNAME}.nwk"
OUTDIR="data/abSENSE_HMM/$DATASETNAME"

OGID="MOG0001047"

mkdir -p $OUTDIR
mkdir -p $OUTDIR/results
mkdir -p $OUTDIR/results_per_species
mkdir -p $OUTDIR/figures

for LOO_HMMSEARCH_FILENAME in $HMM_LOO_RESULTS_PREFIX/results/hmmsearch_out/${OGID}*_LOO_hmmsearch.out
  do
    PHROGID=${LOO_HMMSEARCH_FILENAME%_LOO_hmmsearch.out}
    PHROGID=${PHROGID##*/}
    Rscript run_absense.R $PHROGID $HMM_LOO_RESULTS_PREFIX $ALERAX_DIR $POSTERIOR_CLADES_DIR $OUTDIR $DATASETNAME
  done
done;
