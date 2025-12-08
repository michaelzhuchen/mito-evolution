### Process species tree

library(tidyverse)
library(ape)
library(phytools)
library(TreeTools)

# Read in inferred species tree
inferred_species_tree <- read.tree(here("data/species_phylogeny/maximum_likelihood_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_partitionfile.nex_MFMERGE_rcluster10.best_scheme_constrained.ncbi.tree.manual.changes.v7.1412taxa.contree"))

# Read in prokspp data
bacteria_tax_collapse <- read.table(here("data/taxonomy", "bacteria_taxids_collapse_for_ncbi_taxonomy_tree_taxids_manual.changes.v6.tsv"), sep="\t", header=TRUE)
archaea_tax_collapse <- read.table(here("data/taxonomy", "archaea_taxids_collapse_for_ncbi_taxonomy_tree_taxids_manual.changes.v6.tsv"), sep="\t", header=TRUE)
bacteria_tax_seprows <- bacteria_tax_collapse %>% separate_rows(taxids_for_tree, sep=":1,")
bacteria_tax_seprows$taxids_for_tree <- gsub(":1", "", bacteria_tax_seprows$taxids_for_tree, fixed=TRUE)
archaea_tax_seprows <- archaea_tax_collapse %>% separate_rows(taxids_for_tree, sep=":1,")
archaea_tax_seprows$taxids_for_tree <- gsub(":1", "", archaea_tax_seprows$taxids_for_tree, fixed=TRUE)
bacteria_tax_seprows$kingdom <- gsub("\\/", "", gsub(" ", "", bacteria_tax_seprows$kingdom))
archaea_tax_seprows$superfamily <- gsub("\\/", "", gsub(" ", "", archaea_tax_seprows$superfamily))
colnames(bacteria_tax_seprows)[1] <- "prok_group"
colnames(archaea_tax_seprows)[1] <- "prok_group"
prok_tax_seprows <- rbind(bacteria_tax_seprows, archaea_tax_seprows)

# Rename prokspp tips to prokgroup
inferred_species_tree$tip.label[which(inferred_species_tree$tip.label %in% prok_tax_seprows$taxids_for_tree)] <- prok_tax_seprows$prok_group[match(inferred_species_tree$tip.label[which(inferred_species_tree$tip.label %in% prok_tax_seprows$taxids_for_tree)], prok_tax_seprows$taxids_for_tree)]

# Assign node labels
inferred_species_tree$node.label <- paste0("Node", 1:inferred_species_tree$Nnode)

# Root using phytools:reroot to root at the midpoint of the branch
root_edge_length <- inferred_species_tree$edge.length[inferred_species_tree$edge[,2] == Ntip(inferred_species_tree)+375]
inferred_species_tree_rooted <- reroot(inferred_species_tree, Ntip(inferred_species_tree)+375, position = root_edge_length/2)

# Collapse identical tip labels
tip_labels <- inferred_species_tree_rooted$tip.label
duplicate_labels <- tip_labels[duplicated(tip_labels)]

for (label in unique(duplicate_labels)) {
  print(label)
  # Find the indices of the tips with the same label
  indices <- which(inferred_species_tree_rooted$tip.label == label)
  
  # Compute average length from MRCA to tips
  mrca_index <- get_mrca_of_set(inferred_species_tree_rooted, indices)
  distances_to_mrca <- get_pairwise_distances(inferred_species_tree_rooted, indices, rep(mrca_index, length(indices)), as_edge_counts=FALSE, check_input=TRUE)
  avg_length <- mean(distances_to_mrca)
  print(avg_length)
  
  # Remove the original tips
  inferred_species_tree_rooted <- drop.tip(inferred_species_tree_rooted, indices[-1])
  
  # Update the remaining tip's branch length
  remaining_index <- which(inferred_species_tree_rooted$tip.label == label)
  edge_index <- match(remaining_index, inferred_species_tree_rooted$edge[, 2])
  inferred_species_tree_rooted$edge.length[edge_index] <- avg_length
}

# Sanity check the modified tree
plot(inferred_species_tree_rooted)
inferred_species_tree_rooted

inferred_species_tree_rooted$node.label <- paste0("Node", 1:inferred_species_tree_rooted$Nnode)

# Label nodes
consensus_lineage_labels_per_node <- consensus_taxonomies(inferred_species_tree_rooted, tip_taxonomies = uniprot_proteomes_tax$Lineage[match(inferred_species_tree_rooted$tip.label, uniprot_proteomes_tax$tree_id)], delimiter = "; ")
consensus_lineage_labels_per_node <- gsub(".*; ", "", consensus_lineage_labels_per_node)
inferred_species_tree_rooted$node.label <- paste0("Node", 1:inferred_species_tree_rooted$Nnode, "_", consensus_lineage_labels_per_node)
inferred_species_tree_rooted$node.label[30] <- "Node30_Archaea"
inferred_species_tree_rooted$node.label[31] <- "Node31_Archaea"
inferred_species_tree_rooted$node.label[32] <- "Node32_Archaea"
inferred_species_tree_rooted$node.label[33] <- "Node33_Asgardgroup"
inferred_species_tree_rooted$node.label[428] <- "Node428_Opimoda"
inferred_species_tree_rooted$node.label[429] <- "Node429_Amorphea"
inferred_species_tree_rooted$node.label[430] <- "Node430_Obazoa"
inferred_species_tree_rooted$node.label[35] <- "Node35_Diphoda"
inferred_species_tree_rooted$node.label[36] <- "Node36_Diaphorectickes"
inferred_species_tree_rooted$node.label[37] <- "Node37_CAM_Haptista"
inferred_species_tree_rooted$node.label[38] <- "Node38_CAM"
inferred_species_tree_rooted$node.label[39] <- "Node39_Archaeplastida"

# ## Write out
# write.tree(here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample_v2.contree"))
