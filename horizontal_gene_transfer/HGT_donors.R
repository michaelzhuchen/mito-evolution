### HGT analysis

# Load libraries
library(tidyverse)
library(ape)
library(castor)
library(phytools)


## Select focal taxonomic level
selected_tax_level <- "Node34_Eukaryota" # LECA
# selected_tax_level <- "Node39_Archaeplastida" # Plastid control


## Read in data
# Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
uniprot_proteomes_all_tax <- read.table(here("data/taxonomy", "uniprot_bacteria.downsample.level6_tax_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

# Read in Eukaryota parent PhROGs
parent_progs_long <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(parent_progs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
parent_progs_long <- parent_progs_long %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))

# Get mito Eukaryota_parent PhROGs
parent_progs_long_mito_PhROG_ids <- read.table(here("data/reconstruction", "parent_mito_PhROG_ids_corespecies_deeplocmc.min5spp_Eukaryota.parent.2supergroups.min4spp.filter.txt"))$V1

# Read in HGT PhROGs
combined_hgt_phrog_raw <- read.table(here("data/horizontal_gene_transfer", "posterior_clades_HGT_Node34_Eukaryota_parent_wide.tsv"), sep="\t")
colnames(combined_hgt_phrog_raw) <- c("OG_id","label","clade_index","distance_to_root","reference_protein_ids","count","species_overlap","duplications_rec","nonvertical_protein_ids","n_species","n_reference_proteins","fraction_euk_species","mito_localization_prob","PhROG_id","fraction_primary_OG_for_vertical_proteins","self_clade_protein_ids","sister_clade_protein_ids","cousin_clade_protein_ids","grandma_clade_protein_ids","self_clade_prok_species_in_clades","self_clade_n_prok_species","self_clade_n_euk_species","self_clade_fraction_euk_species","sister_clade_prok_species_in_clades","sister_clade_n_prok_species","sister_clade_n_euk_species","sister_clade_fraction_euk_species","cousin_clade_prok_species_in_clades","cousin_clade_n_prok_species","cousin_clade_n_euk_species","cousin_clade_fraction_euk_species","grandma_clade_prok_species_in_clades","grandma_clade_n_prok_species","grandma_clade_n_euk_species","grandma_clade_fraction_euk_species","HGT_self", "HGT_sister","HGT_cousin")

# Convert NAs to empty string to avoid errors
combined_hgt_phrog_raw$self_clade_protein_ids[is.na(combined_hgt_phrog_raw$self_clade_protein_ids)] <- ""
combined_hgt_phrog_raw$sister_clade_protein_ids[is.na(combined_hgt_phrog_raw$sister_clade_protein_ids)] <- ""
combined_hgt_phrog_raw$cousin_clade_protein_ids[is.na(combined_hgt_phrog_raw$cousin_clade_protein_ids)] <- ""
combined_hgt_phrog_raw$grandma_clade_protein_ids[is.na(combined_hgt_phrog_raw$grandma_clade_protein_ids)] <- ""

# Convert infinite values to 0
combined_hgt_phrog_raw$sister_clade_fraction_euk_species[is.infinite(combined_hgt_phrog_raw$sister_clade_fraction_euk_species)] <- 0
combined_hgt_phrog_raw$cousin_clade_fraction_euk_species[is.infinite(combined_hgt_phrog_raw$cousin_clade_fraction_euk_species)] <- 0

# Function to find largest parenthetical element in comma-separated string
get_largest_elems <- function(input_string, min_threshold) {
  if (is.na(input_string)) {
    return("")
  }
  if (input_string == "") {
    return("")
  }
  
  # Split the string into individual elements
  elements <- strsplit(input_string, ",")[[1]]
  
  # Remove empty elements
  elements <- elements[elements != ""]
  
  # Extract the names and the numeric values inside parentheses
  name_values <- gsub(".*\\(([^)]+)\\)", "\\1", elements)  # Extract numbers
  names <- gsub("\\(.*\\)", "", elements)  # Extract names
  
  # Convert the numeric values to a numeric vector
  num_values <- as.numeric(name_values)
  
  # Find the maximum number
  max_value <- max(num_values)
  
  # Extract the names corresponding to the maximum value
  max_names <- names[num_values == max_value]
  
  if (max_value >= min_threshold) {
    return(paste0(paste0(max_names, "(", max_value, ")"), collapse=","))
  } else {
    return("")
  }
}


## Identify HGT donors

# Set minimum relative fraction species coverage per prokaryotic group. Require enrichment equal to at least the baseline frequency
fraction_species_coverage_threshold <- 1

combined_hgt_phrog <- combined_hgt_phrog_raw

## Filter for LECA clades with prokaryotic origin
origin_table <- read.table(here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t")
colnames(origin_table) <- c("OG_id", "n_euk_species_largest_euk_clade", "n_prok_species_largest_prok_clade", "LCA_node", "LCA_node_euks", "origin_domain", "dropped_tips")
origin_table$OG_id <- gsub(".faa_clipkit.gappy.msa", "", origin_table$OG_id, fixed=TRUE)
prok_origin_OG_ids <- origin_table$OG_id[origin_table$origin_domain == "Prokaryote"]

combined_hgt_phrog_prok_origin <- combined_hgt_phrog %>% filter(OG_id %in% prok_origin_OG_ids)

# Get PhROGs at/above the selected taxonomic level
species_tree <- read.tree(here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample_v2.contree"))
species_tree_subtree <- get_subtree_at_node(species_tree, selected_tax_level)$subtree
species_tree_euks_subtree <- get_subtree_at_node(species_tree, "Node34_Eukaryota")$subtree
species_tree_labels <- c(species_tree$tip.label, species_tree$node.label)
species_tree_euks_labels <- c(species_tree_euks_subtree$tip.label, species_tree_euks_subtree$node.label)
species_tree_subtree_labels <- c(species_tree_subtree$tip.label, species_tree_subtree$node.label)

# Get ancestors on path to root
curr_tax_label <- selected_tax_level
curr_edge_index <- which(species_tree_labels == curr_tax_label)
ancestral_tax_labels <- c(curr_tax_label)
while (curr_edge_index != find_root(species_tree)) {
  curr_edge_index <- species_tree$edge[species_tree$edge[,2] == curr_edge_index,1]
  curr_tax_label <- species_tree_labels[curr_edge_index]
  ancestral_tax_labels <- c(ancestral_tax_labels, curr_tax_label)
  curr_edge_index <- which(species_tree_labels == curr_tax_label)
}
ancestral_tax_labels <- ancestral_tax_labels[which(!ancestral_tax_labels %in% species_tree_subtree_labels)]
ancestral_tax_labels <- c(ancestral_tax_labels, selected_tax_level)

parent_progs_long_selected_tax_level <- parent_progs_long %>% filter(taxid %in% species_tree_subtree$tip.label)
parent_progs_long_selected_tax_level_PhROG_ids <- unique(parent_progs_long_selected_tax_level$PROG_id)
combined_hgt_phrog_ancestor_prok_origin <- combined_hgt_phrog_prok_origin %>% rowwise() %>% filter(PhROG_id %in% parent_progs_long_selected_tax_level_PhROG_ids) %>% filter(label %in% ancestral_tax_labels)

## Compute fraction of species per group, (adjusted for number of species in group) present in sister clade
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin %>% rowwise() %>% mutate(self_sister_clade_prok_species_in_clades = sister_clade_prok_species_in_clades)

combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% separate_rows(self_sister_clade_prok_species_in_clades, sep=",") %>% filter(self_sister_clade_prok_species_in_clades != "") %>% group_by(PhROG_id) %>% filter(!duplicated(self_sister_clade_prok_species_in_clades))
combined_hgt_phrog_ancestor_prok_origin_species_counts$tree_id <- uniprot_proteomes_all_tax$tree_id[match(combined_hgt_phrog_ancestor_prok_origin_species_counts$self_sister_clade_prok_species_in_clades, uniprot_proteomes_all_tax$TaxId)]
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% group_by(PhROG_id, tree_id) %>% summarize(n_species = n())

uniprot_proteomes_all_tax_prokaryote_species_counts <- uniprot_proteomes_all_tax %>% filter(domain %in% c("Bacteria", "Archaea")) %>% group_by(tree_id) %>% summarize(total_species_count = n())
combined_hgt_phrog_ancestor_prok_origin_species_counts$total_species_count <- uniprot_proteomes_all_tax_prokaryote_species_counts$total_species_count[match(combined_hgt_phrog_ancestor_prok_origin_species_counts$tree_id, uniprot_proteomes_all_tax_prokaryote_species_counts$tree_id)]
# Fraction of clade corrected for prokgroup size
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% group_by(PhROG_id) %>% mutate(raw_fraction_species_coverage = (n_species / sum(n_species)), fraction_group_species_of_total_species = total_species_count / sum(uniprot_proteomes_all_tax_prokaryote_species_counts$total_species_count)) %>% mutate(fraction_species_coverage = raw_fraction_species_coverage / fraction_group_species_of_total_species)

combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts[order(combined_hgt_phrog_ancestor_prok_origin_species_counts$fraction_species_coverage, decreasing=TRUE),]
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% rowwise() %>% mutate(combined_label = paste0(tree_id, "(", fraction_species_coverage, ")"))
combined_hgt_phrog_ancestor_prok_origin_species_counts_summary <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% group_by(PhROG_id) %>% summarize(HGT_self_sister = paste0(combined_label, collapse=","))
combined_hgt_phrog_ancestor_prok_origin$HGT_self_sister <- combined_hgt_phrog_ancestor_prok_origin_species_counts_summary$HGT_self_sister[match(combined_hgt_phrog_ancestor_prok_origin$PhROG_id, combined_hgt_phrog_ancestor_prok_origin_species_counts_summary$PhROG_id)]
combined_hgt_phrog_ancestor_prok_origin <- combined_hgt_phrog_ancestor_prok_origin %>% filter(!is.na(HGT_self_sister)) %>% filter(HGT_self_sister != "")

# Assign most likely prokgroup by highest enrichment in sister clade
combined_hgt_phrog_ancestor_prok_origin <- combined_hgt_phrog_ancestor_prok_origin %>% rowwise() %>% mutate(HGT_self_sister = get_largest_elems(HGT_self_sister, fraction_species_coverage_threshold))

# Filter by absence of eukaryotes in sister+cousin clade as evidence of recent HGT from prokaryotes
combined_hgt_phrog_ancestor_prok_origin <- combined_hgt_phrog_ancestor_prok_origin %>% filter(sister_clade_fraction_euk_species == 0 & cousin_clade_fraction_euk_species == 0)

# Label with mito/nonmito localization
combined_hgt_phrog_ancestor_prok_origin <- combined_hgt_phrog_ancestor_prok_origin %>% mutate(bool_in_mito = PhROG_id %in% parent_progs_long_mito_PhROG_ids)

# Convert to long format
combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long <- combined_hgt_phrog_ancestor_prok_origin %>% separate_rows(HGT_self_sister, sep=",") %>% rowwise() %>% mutate(taxon = gsub("\\(.*\\)", "", HGT_self_sister), fraction_species_coverage = as.numeric(gsub(".*\\(([^)]+)\\)", "\\1", HGT_self_sister))) %>% filter(fraction_species_coverage >= fraction_species_coverage_threshold)

# Subset and rename columns
leca_HGT_long <- leca_HGT_long[,c("OG_id", "PhROG_id", "count", "fraction_euk_species", "sister_clade_fraction_euk_species", "cousin_clade_fraction_euk_species", "HGT_sister", "HGT_self_sister", "taxon", "bool_in_mito")]
colnames(leca_HGT_long) <- c("OG_ID", "PhROG_ID", "Support (conditional_clade_probability x 100)", "Fraction_eukaryotic_species_in_clade", "Fraction_eukaryotic_species_in_sister_clade", "Fraction_eukaryotic_species_in_cousin_clade", "Prokaryotic_groups_in_sister_clade (fraction_of_prokaryotic_group_species)", "Top_prokaryotic_group_in_sister_clade (enrichment_score)", "Top_prokaryotic_group_in_sister_clade", "In_reconstructed_ancestral_mitochondria")

## Write out
# write.table(combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long, here("data/horizontal_gene_transfer", paste0("HGT_", selected_tax_level, "_long.tsv")), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


