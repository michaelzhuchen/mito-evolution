#!/bin/sh


families_file_template="families_template_pruned.ufboot.rooted.txt"
species_tree="data/species_phylogeny/processed_species_tree/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample.contree"
base_output_directory="output_directory"
parametrization_file="alerax_model_parameterization_Node34_Eukaryota.txt"
fraction_missing_file="fraction_missing_euk203spp_prokgroups_alerax.txt"


for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
	MSAFILE="${OGID}.faa_clipkit.gappy.msa"

	OUTPUT_DIR="${base_output_directory}/${OGID}"

	mkdir -p $OUTPUT_DIR

	# Create families file for this OG
	families_file="${OUTPUT_DIR}/families_${OGID}.txt"
	sed "s|POSTERIORTREESDIR|${BASE_DIR}|g" "$families_file_template" > "$families_file"
	sed -i "s/MSAFILE/${MSAFILE}/g" "$families_file"

	# Eukaryote-specific rates per family and rooted gene trees.
	alerax -f $families_file -s $species_tree -p $OUTPUT_DIR --rec-opt GRADIENT --gene-tree-samples 100 --model-parametrization $parametrization_file --fraction-missing-file $fraction_missing_file --seed 42 --gene-tree-rooting ROOTED
done;
