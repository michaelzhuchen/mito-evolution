#!/bin/bash

N_CORES="16"

DIR="zenodo_directory/species_phylogeny"
MSA_FILE="$DIR/concatenated_msa/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa"

## 1a. Infer optimal partition models without building tree, using rcluster heuristic for speedup
PARTITIONFILE="$DIR/concatenated_msa/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_partitionfile.nex"
OUTPUT_PREFIX="$DIR/maximum_likelihood_species_tree/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_partitionfile.nex_MFMERGE_rcluster10"
/home/cmichael/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $MSA_FILE -m MF+MERGE -cmin 4 -cmax 12 -p $PARTITIONFILE -nt $N_CORES -safe -seed 42 -rcluster 10 --prefix $OUTPUT_PREFIX

## 1b. Run with partition model with pre-existing model selection
PARTITIONFILE="$DIR/maximum_likelihood_species_tree/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_partitionfile.nex_MFMERGE_rcluster10.best_scheme.nex"
OUTPUT_PREFIX="$DIR/maximum_likelihood_species_tree/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_partitionfile.nex_MFMERGE_rcluster10.best_scheme_constrained.ncbi.tree.manual.changes.v7.1412taxa"
CONSTRAINT_TREE="$DIR/guide_tree/ncbi_taxonomy_tree_taxids_manual.changes.v7_prokspp_no.node.labels_collapse.single.nodes_1412taxa.nwk"
/home/cmichael/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $MSA_FILE -p $PARTITIONFILE -g $CONSTRAINT_TREE --prefix $OUTPUT_PREFIX -nt $N_CORES -bb 1000 -alrt 1000 -safe -seed 42 -bnni
