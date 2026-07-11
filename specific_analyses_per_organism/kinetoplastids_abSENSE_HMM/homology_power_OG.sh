#!/bin/sh

DATASETNAME="species_tree_1"
ALERAX_DIR="reconciled_trees_species.tree.1"
POSTERIOR_CLADES_DIR="reconciled_trees_posterior_clades_species.tree.1"
HMM_LOO_RESULTS_PREFIX="data/abSENSE_HMM/$DATASETNAME"
SPECIES_TREE_FILENAME="data/species_phylogeny/processed_species_tree/${DATASETNAME}.nwk"
OUTDIR="data/abSENSE_HMM/$DATASETNAME"

OGID="MOG0001047"

Rscript abSENSE_HMM/run_absense.R $OGID $HMM_LOO_RESULTS_PREFIX $ALERAX_DIR $POSTERIOR_CLADES_DIR $OUTDIR $DATASETNAME
