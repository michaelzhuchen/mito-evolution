# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(castor))

# Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
base_dir <- args[1]
out_dir <- args[2]
OG_id <- args[3]
selected_tax_levels_file <- args[4]
BOOL_verbose <- FALSE

# Focal taxonomic levels for which to generate PhROGs
selected_tax_levels <- read.table(here("phylogenetically_resolved_orthogroups", "selected_tax_levels_for_phrogs_wholeproteome.txt"), header=FALSE)$V1

filename <- paste0(OG_id, "_euk_monophyletic_clades.tsv")

clade_support_threshold <- 50
species_overlap_threshold <- 0.1
species_overlap_fraction_species_coverage_threshold <- 0.1
clade_purity_threshold <- 0.9
species_overlap_support_threshold <- 0.5
duplications_rec_support_threshold <- 0
n_euk_species_in_dataset <- 203
n_archaea_species_in_dataset <- 337
n_bacteria_species_in_dataset <- 1737
fraction_euk_species_threshold <- round(2 * (n_euk_species_in_dataset) / (n_euk_species_in_dataset + n_archaea_species_in_dataset + n_bacteria_species_in_dataset), digits=2)

# euk203spp_prokgroups
species_tree <- read.tree(here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample_v2.contree"))
species_tree_euks_subtree <- get_subtree_at_node(species_tree, "Node34_Eukaryota")$subtree
species_tree_labels <- c(species_tree$tip.label, species_tree$node.label)
species_tree_euks_labels <- c(species_tree_euks_subtree$tip.label, species_tree_euks_subtree$node.label)
# Calculate node depths
depths_df <- data.frame(label = c(species_tree$tip.label, species_tree$node.label), distance=get_all_distances_to_root(species_tree, as_edge_count=TRUE))

prok_euk_parent_labels <- c("Node1_cellular_organisms", "Node30_Archaea", "Node31_Archaea", "Node32_Archaea", "Node33_Asgardgroup")

# Read in fusion proteins primary OG assignments to designate primary OGs
fusion_protein_primary_disjoint_OG_mapping <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "euk673spp_prokgroups_OG_all_merged.hmm.foldseek_add.species.singleton.mtDNA_50AA_expect1e-10_fusion.protein_to_OG.txt"), sep="\t", header=TRUE)

## Preprocess posterior clades
if (BOOL_verbose) {
  print("Preprocessing clades...")
}
# Read in posterior clades 
monophyletic_clades <- read.table(paste0(base_dir, "/", filename), sep="\t")
colnames(monophyletic_clades) <- c("OG_id", "label", "distance_to_root", "reference_protein_ids", "nonvertical_protein_ids", "count", "species_overlap", "species_overlap_support", "species_overlap_taxids", "duplications_rec", "mito_localization_prob", "n_species", "n_reference_proteins") # vtesting

# euk203spp_prokgroups: Map species tree v1 nodelabels (manual assignments were incorrect due to off by one) to v2 nodelabels
monophyletic_clades$label[which(monophyletic_clades$label == "Node138_Opimoda")] <- "Node138_Naegleria"
monophyletic_clades$label[which(monophyletic_clades$label == "Node139_Amorphea_CRuMs")] <- "Node139_Opimoda"
monophyletic_clades$label[which(monophyletic_clades$label == "Node140_Amorphea")] <- "Node140_Amorphea_CRuMs"
monophyletic_clades$label[which(monophyletic_clades$label == "Node141_Obazoa")] <- "Node141_Amorphea"
monophyletic_clades$label[which(monophyletic_clades$label == "Node142_Eukaryota")] <- "Node142_Obazoa"

# Map lta2019 to ltaref
# Read in Lta mapping for Lta2019 to Ltaref
lta_mapping <- read.delim(here("data/orthogroups/idmapping", "map.ltaref.lta.exact.txt"), header=FALSE)
colnames(lta_mapping) <- c("ltaref_id", "lta2019_id")
lta_mapping$ltaref_id <- paste0("5689_", lta_mapping$ltaref_id)
lta_mapping$lta2019_id[lta_mapping$lta2019_id == ""] <- NA
lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)] <- paste0("5689_", lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])

## Map new Lta accessions and keep only the identical proteins in ltaref (5689)
# Keep all Lta proteins that are identical in ltaref
filter_items <- function(item_string, allowed_items) {
  if (is.na(item_string)) {
    return("")
  }
  items <- unlist(strsplit(item_string, ","))
  items_other <- items[!grepl("^5689_GET", items)]
  items_selected <- items[grepl("^5689_GET", items)]
  filtered_items <- intersect(items_selected, allowed_items)
  filtered_items <- c(filtered_items, items_other)
  paste(filtered_items, collapse = ",")
}
# Replace Lta2019 protein ids with ltaref protein ids
replace_items <- function(item_string, mapping_table) {
  if (is.na(item_string)) {
    return("")
  }
  items <- unlist(strsplit(item_string, ","))
  items_other <- items[!grepl("^5689_GET", items)]
  items_selected <- items[grepl("^5689_GET", items)]
  mapped_items <- mapping_table$ltaref_id[which(mapping_table$lta2019_id %in% items_selected)]
  mapped_items <- c(mapped_items, items_other)
  paste(mapped_items, collapse = ",")
}

monophyletic_clades$reference_protein_ids <- sapply(monophyletic_clades$reference_protein_ids, filter_items, lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
monophyletic_clades$reference_protein_ids <- sapply(monophyletic_clades$reference_protein_ids, replace_items, lta_mapping)
monophyletic_clades$nonvertical_protein_ids <- sapply(monophyletic_clades$nonvertical_protein_ids, filter_items, lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
monophyletic_clades$nonvertical_protein_ids <- sapply(monophyletic_clades$nonvertical_protein_ids, replace_items, lta_mapping)

# Filter by clade support
monophyletic_clades_filter <- monophyletic_clades %>% filter(count >= clade_support_threshold)

# Refine taxonomic labels for mixed prok+euk clades
# Read in the prokaryote mmseqs2 cluster data
prokaryote_mmseqs2_clusters_taxids <- read.table(here("data/downsample_prokaryotes", "prokaryote_mmseqs2_clusters_taxids.tsv"), sep="\t", quote="", header=TRUE)
all_proteins_remove_first_underscore <- sub("^[^_]*_", "", unlist(strsplit(monophyletic_clades_filter$reference_protein_ids, split=",")))
prokaryote_mmseqs2_clusters_taxids <- prokaryote_mmseqs2_clusters_taxids[which(prokaryote_mmseqs2_clusters_taxids$rep_seq %in% all_proteins_remove_first_underscore),]
# Get prokaryote proteins and species counts
prok_proteins_in_clades_raw <- lapply(strsplit(monophyletic_clades_filter$reference_protein_ids, split=","), function(lst) lst[!gsub("_.*", "", lst) %in% species_tree_euks_subtree$tip.label])
prok_proteins_in_clades <- lapply(prok_proteins_in_clades_raw, function(lst) sub("^[^_]*_", "", lst))
prok_species_in_clades <- unlist(lapply(prok_proteins_in_clades, function(lst) {
  paste0(sort(unique(unlist(strsplit(prokaryote_mmseqs2_clusters_taxids$unique_member_seq_taxids[which(prokaryote_mmseqs2_clusters_taxids$rep_seq %in% lst)], split=",")))), collapse=",")
}))
monophyletic_clades_filter$n_prok_species <- unlist(lapply(prok_species_in_clades, function(lst) length(which(!is.na(unlist(strsplit(lst, split=",")))))))
monophyletic_clades_filter <- monophyletic_clades_filter %>% rowwise() %>% mutate(n_euk_species = sum(species_tree_euks_subtree$tip.label %in% gsub("_.*", "", unlist(strsplit(reference_protein_ids, split=",")))))
monophyletic_clades_filter <- monophyletic_clades_filter %>% rowwise() %>% mutate(fraction_euk_species = n_euk_species / (n_euk_species + n_prok_species))

# Rename majority euk clades to euk labels
relabel_euk_indexes <- which(monophyletic_clades_filter$fraction_euk_species > clade_purity_threshold & monophyletic_clades_filter$label %in% prok_euk_parent_labels)
if (length(relabel_euk_indexes) > 0) {
  for (relabel_euk_index in relabel_euk_indexes) {
    curr_reference_species_ids <- unique(gsub("_.*", "", unlist(strsplit(monophyletic_clades_filter$reference_protein_ids[relabel_euk_index], split=","))))
    curr_reference_species_ids <- curr_reference_species_ids[curr_reference_species_ids %in% species_tree_euks_subtree$tip.label]
    lca_index <- get_mrca_of_set(species_tree, curr_reference_species_ids)
    if (lca_index <= Ntip(species_tree)) {
      curr_species_tree_label <- species_tree$tip.label[lca_index]
    } else {
      curr_species_tree_label <- species_tree$node.label[lca_index - Ntip(species_tree)]
    }
    if (BOOL_verbose) {
      print(paste0("Original label: ", monophyletic_clades_filter$label[relabel_euk_index]))
      print(paste0("New label: ", curr_species_tree_label))
    }
    monophyletic_clades_filter$label[relabel_euk_index] <- curr_species_tree_label
    monophyletic_clades_filter$distance_to_root[relabel_euk_index] <- depths_df$distance[match(monophyletic_clades_filter$label[relabel_euk_index], depths_df$label)]
  }
}

if (BOOL_verbose) {
  write.table(monophyletic_clades_filter, paste0(out_dir, "/", OG_id, "_relabeled_clades.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
}


## Filter for eukaryote-enriched clades to identify independent gains into eukaryotes
monophyletic_clades_filter_selected_tax_level <- monophyletic_clades_filter %>% filter(fraction_euk_species > fraction_euk_species_threshold)

# Generate PhROGs at each taxonomic level
for (selected_tax_level in selected_tax_levels) {
  if (BOOL_verbose) {
    print(selected_tax_level)
  }
  original_selected_tax_level <- selected_tax_level
  if (grepl("_parent$", selected_tax_level)) {
    BOOL_split_duplications_at_selected_tax_level <- FALSE
    selected_tax_level <- gsub("_parent$", "", selected_tax_level)
  } else {
    BOOL_split_duplications_at_selected_tax_level <- TRUE
  }
  
  curr_out_dir <- paste0(out_dir, "/", original_selected_tax_level)
  dir.create(curr_out_dir, showWarnings = FALSE)
  
  # Check for completion
  if (file.exists(paste0(curr_out_dir, "/", OG_id, "_", original_selected_tax_level, "_PhROGs_long.tsv"))) {
    if (BOOL_verbose) {
      print("Already completed")
    }
    next
  }
  
  species_tree_subtree <- get_subtree_at_node(species_tree, selected_tax_level)$subtree
  species_tree_subtree_labels <- c(species_tree_subtree$tip.label, species_tree_subtree$node.label)
  # Split or lump duplications that occur in the selected tax level. Set to TRUE if want to get genes in the LAST common ancestor
  if (BOOL_split_duplications_at_selected_tax_level) {
    species_tree_subtree_labels_selected <- species_tree_subtree_labels[which(species_tree_subtree_labels != selected_tax_level)]
  } else {
    species_tree_subtree_labels_selected <- species_tree_subtree_labels
  }
  
  # Get daughter subtrees
  node_index <- which(species_tree_subtree_labels == selected_tax_level)
  daughter_indexes <- species_tree_subtree$edge[which(species_tree_subtree$edge[,1] == node_index),2]
  daughter_node_labels <- c()
  subtrees <- c()
  subtree_tip_labels <- c()
  for (daughter_index in daughter_indexes) {
    if (daughter_index > Ntip(species_tree_subtree)) {
      # If it's a node
      daughter_node_label <- species_tree_subtree$node.label[daughter_index - Ntip(species_tree_subtree)]
      daughter_node_labels <- c(daughter_node_labels, daughter_node_label)
      subtree_tip_labels <- c(subtree_tip_labels, list(get_subtree_at_node(species_tree_subtree, daughter_index - Ntip(species_tree_subtree))$subtree$tip.label))
      if (length(subtrees) == 0) {
        subtrees <- c(get_subtree_at_node(species_tree_subtree, daughter_index - Ntip(species_tree_subtree))$subtree)
      } else {
        subtrees <- c(subtrees, get_subtree_at_node(species_tree_subtree, daughter_index - Ntip(species_tree_subtree))$subtree)
      }
    } else {
      # If it's a tip, extract the taxid
      daughter_node_label <- gsub("_.*", "", species_tree_subtree$tip.label[daughter_index])
      daughter_node_labels <- c(daughter_node_labels, daughter_node_label)
      if (length(subtrees) == 0) {
        subtrees <- c(keep.tip(species_tree_subtree, species_tree_subtree$tip.label[daughter_index]))
      } else {
        subtrees <- c(subtrees, keep.tip(species_tree_subtree, species_tree_subtree$tip.label[daughter_index]))
      }
      subtree_tip_labels <- c(subtree_tip_labels, list(species_tree_subtree$tip.label[daughter_index]))
    }
  }
  
  # Find largest monophyletic clades per protein
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level
  n_clades_initial <- nrow(monophyletic_clades_filter_selected_tax_level)
  monophyletic_clades_filter_selected_tax_level_sep_rows$clade_index <- 1:n_clades_initial
  
  # Keep just the species beneath the selected tax level, but preserve rows not in tax level for finding duplications later.
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows[order(monophyletic_clades_filter_selected_tax_level_sep_rows$n_reference_proteins, 1/monophyletic_clades_filter_selected_tax_level_sep_rows$distance_to_root, decreasing=TRUE),]
  monophyletic_clades_filter_selected_tax_level_sep_rows$reference_protein_ids_total <- monophyletic_clades_filter_selected_tax_level_sep_rows$reference_protein_ids
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% separate_rows(reference_protein_ids, sep=",")
  monophyletic_clades_filter_selected_tax_level_sep_rows_intaxlevel <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(reference_protein_ids != "") %>% filter(!is.na(reference_protein_ids)) %>% filter(gsub("_.*", "", reference_protein_ids) %in% species_tree_subtree$tip.label)
  # Preserve rows not in tax level for finding duplications later. keep just one per original row.
  monophyletic_clades_filter_selected_tax_level_sep_rows_notintaxlevel <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(reference_protein_ids != "") %>% filter(!is.na(reference_protein_ids)) %>% filter(!gsub("_.*", "", reference_protein_ids) %in% species_tree_subtree$tip.label) %>% filter(!duplicated(reference_protein_ids_total))
  monophyletic_clades_filter_selected_tax_level_sep_rows <- rbind(monophyletic_clades_filter_selected_tax_level_sep_rows_intaxlevel, monophyletic_clades_filter_selected_tax_level_sep_rows_notintaxlevel)
  
  # Replace NA entries with empty string to avoid errors in summarize
  monophyletic_clades_filter_selected_tax_level_sep_rows$nonvertical_protein_ids[is.na(monophyletic_clades_filter_selected_tax_level_sep_rows$nonvertical_protein_ids)] <- ""
  monophyletic_clades_filter_selected_tax_level_sep_rows$species_overlap_taxids[is.na(monophyletic_clades_filter_selected_tax_level_sep_rows$species_overlap_taxids)] <- ""
  
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% group_by(OG_id, label, clade_index, reference_protein_ids_total) %>% summarize(distance_to_root = mean(distance_to_root), reference_protein_ids = paste0(sort(unique(unlist(strsplit(reference_protein_ids, split=",")))), collapse=","), count = mean(count), species_overlap = mean(species_overlap), species_overlap_support = mean(species_overlap_support), duplications_rec = mean(duplications_rec), species_overlap_taxids = paste0(sort(unique(unlist(strsplit(species_overlap_taxids, split=",")))), collapse=","), nonvertical_protein_ids = paste0(sort(unique(unlist(strsplit(nonvertical_protein_ids, split=",")))), collapse=","), n_species_total = mean(n_species), n_reference_proteins_total = mean(n_reference_proteins), n_species = length(unique(gsub("_.*", "", unlist(strsplit(reference_protein_ids, split=","))))), n_reference_proteins = length(unlist(strsplit(reference_protein_ids, split=","))), n_prok_species = mean(n_prok_species), n_euk_species = mean(n_euk_species), fraction_euk_species = mean(fraction_euk_species), mito_localization_prob = mean(mito_localization_prob), .groups = "drop")
  
  # If no species under taxonomic level are in the OG, write empty file and skip
  if (nrow(monophyletic_clades_filter_selected_tax_level_sep_rows_intaxlevel) == 0) {
    file.create(paste0(curr_out_dir, "/", OG_id, "_", original_selected_tax_level, "_PhROGs_long.tsv"))
    next
  }
  
  ### Split at duplications upstream or in the current taxonomic level
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% rowwise() %>% mutate(n_species_overlap_in_tax_level = sum(unique(unlist(strsplit(species_overlap_taxids, split=","))) %in% species_tree_subtree$tip.label), taxid_species_overlap_in_tax_level_subtree1 = paste0(unique(unlist(strsplit(species_overlap_taxids, split=",")))[unique(unlist(strsplit(species_overlap_taxids, split=","))) %in% subtree_tip_labels[[1]]], collapse=","), n_species_overlap_in_tax_level_subtree1 = sum(unique(unlist(strsplit(species_overlap_taxids, split=","))) %in% subtree_tip_labels[[1]]), taxid_species_overlap_in_tax_level_subtree2 = paste0(unique(unlist(strsplit(species_overlap_taxids, split=",")))[unique(unlist(strsplit(species_overlap_taxids, split=","))) %in% subtree_tip_labels[[2]]], collapse=","), n_species_overlap_in_tax_level_subtree2 = sum(unique(unlist(strsplit(species_overlap_taxids, split=","))) %in% subtree_tip_labels[[2]]))
  
  # Infer duplication timing labels
  monophyletic_clades_filter_selected_tax_level_sep_rows$duplication_label <- monophyletic_clades_filter_selected_tax_level_sep_rows$label
  
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
  ancestral_tax_labels <- ancestral_tax_labels[which(!ancestral_tax_labels %in% species_tree_subtree_labels_selected)]
  
  monophyletic_clades_filter_post_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(duplication_label %in% species_tree_subtree_labels_selected)
  monophyletic_clades_filter_preceding_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(duplication_label %in% ancestral_tax_labels)
  monophyletic_clades_filter_preceding_tax_level_sep_rows <- monophyletic_clades_filter_preceding_tax_level_sep_rows[order(monophyletic_clades_filter_preceding_tax_level_sep_rows$n_reference_proteins_total, monophyletic_clades_filter_preceding_tax_level_sep_rows$n_species_total, decreasing=TRUE),]
  
  monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications <- monophyletic_clades_filter_preceding_tax_level_sep_rows %>% filter(species_overlap >= species_overlap_threshold & species_overlap_support >= species_overlap_support_threshold & duplications_rec >= duplications_rec_support_threshold)
  monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications <- monophyletic_clades_filter_preceding_tax_level_sep_rows %>% filter(!(species_overlap >= species_overlap_threshold & species_overlap_support >= species_overlap_support_threshold & duplications_rec >= duplications_rec_support_threshold))
  
  if (BOOL_verbose) {
    write.table(monophyletic_clades_filter_preceding_tax_level_sep_rows, paste0(curr_out_dir, "/", OG_id, "_clades.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    write.table(monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications, paste0(curr_out_dir, "/", OG_id, "_clades_preceding_duplications.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    write.table(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications, paste0(curr_out_dir, "/", OG_id, "_clades_preceding_noduplications.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
  }
  
  # Remove ancestral rows that precede duplications
  rows_with_descending_duplications <- c()
  if (nrow(monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications) > 0 & nrow(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications) > 0) {
    for (j in 1:nrow(monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications)) {
      # Get all protein ids from the duplicated row, including vertical and nonvertical. Important to consider all proteins for finding the parent rows
      curr_duplication_protein_ids_vertical <- unlist(strsplit(monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications$reference_protein_ids_total[j], split=","))
      curr_duplication_protein_ids_nonvertical <- unlist(strsplit(monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications$nonvertical_protein_ids[j], split=","))
      curr_duplication_protein_ids <- c(curr_duplication_protein_ids_vertical, curr_duplication_protein_ids_nonvertical)
      curr_duplication_n_species_total <- monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications$n_species_total[j]
      curr_duplication_species_overlap <- monophyletic_clades_filter_preceding_tax_level_sep_rows_duplications$species_overlap[j]
      curr_duplication_n_species_duplicated <- curr_duplication_species_overlap * curr_duplication_n_species_total
      
      rows_with_descending_duplications_candidates <- which(unlist(lapply(paste(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications$reference_protein_ids_total, monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications$nonvertical_protein_ids, sep=","), function(lst) {
        protein_ids <- unlist(strsplit(lst, split=","))
        
        # Require that all duplicated proteins are all found in ancestral row
        sum(curr_duplication_protein_ids %in% protein_ids) == length(curr_duplication_protein_ids)
      })))
      
      # Filter by fraction of species overlapping over number of species in row
      if (length(rows_with_descending_duplications_candidates) > 0) {
        monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates <- monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications[rows_with_descending_duplications_candidates,]
        monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates$index <- rows_with_descending_duplications_candidates
        monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates <- monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates %>% rowwise() %>% mutate(fraction_species_overlap = curr_duplication_n_species_duplicated / n_species_total)
        monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates <- monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates %>% filter(fraction_species_overlap > species_overlap_fraction_species_coverage_threshold)
        rows_with_descending_duplications_candidates <- monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_candidates$index
      }
      rows_with_descending_duplications <- c(rows_with_descending_duplications, rows_with_descending_duplications_candidates)
    }
    
    rows_with_descending_duplications <- unique(rows_with_descending_duplications)
    rows_without_descending_duplications <- 1:nrow(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications)
    rows_without_descending_duplications <- rows_without_descending_duplications[!rows_without_descending_duplications %in% rows_with_descending_duplications]
    monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_noparent <- monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications[rows_without_descending_duplications,]
    monophyletic_clades_filter_selected_tax_level_sep_rows <- rbind(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications_noparent, monophyletic_clades_filter_post_tax_level_sep_rows)
  } else if (nrow(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications) > 0)  {
    monophyletic_clades_filter_selected_tax_level_sep_rows <- rbind(monophyletic_clades_filter_preceding_tax_level_sep_rows_noduplications, monophyletic_clades_filter_post_tax_level_sep_rows)
  } else {
    monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_post_tax_level_sep_rows
  }
  
  
  ## Get largest clades for each protein
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows[order(monophyletic_clades_filter_selected_tax_level_sep_rows$n_reference_proteins, 1/monophyletic_clades_filter_selected_tax_level_sep_rows$distance_to_root, decreasing=TRUE),]
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% separate_rows(reference_protein_ids, sep=",")
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(reference_protein_ids != "") %>% filter(!is.na(reference_protein_ids))
  # If no rows remaining, skip
  if (nrow(monophyletic_clades_filter_selected_tax_level_sep_rows) == 0) {
    next
  }
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(gsub("_.*", "", unlist(strsplit(reference_protein_ids, split=","))) %in% species_tree_subtree$tip.label)
  monophyletic_clades_filter_selected_tax_level_sep_rows_largest <- monophyletic_clades_filter_selected_tax_level_sep_rows[!duplicated(monophyletic_clades_filter_selected_tax_level_sep_rows$reference_protein_ids),]
  
  # Find the most specific clades that contain all the species within the selected tax level. i.e. minimum subclade that spans all species within the selected tax level.
  retain_clade_indexes <- c()
  for (curr_clade_index in unique(monophyletic_clades_filter_selected_tax_level_sep_rows_largest$clade_index)) {
    # Get selected tax level species protein ids for largest clade (including nonvertical)
    monophyletic_clades_filter_selected_tax_level_sep_rows_curr <- monophyletic_clades_filter_selected_tax_level_sep_rows_largest %>% filter(clade_index == curr_clade_index)
    unique_nonvertical_protein_ids <- unique(unlist(strsplit(monophyletic_clades_filter_selected_tax_level_sep_rows_curr$nonvertical_protein_ids, split=",")))
    unique_nonvertical_protein_ids <- unique_nonvertical_protein_ids[gsub("_.*", "", unique_nonvertical_protein_ids) %in% species_tree_subtree$tip.label]
    unique_reference_protein_ids <- unique(monophyletic_clades_filter_selected_tax_level_sep_rows_curr$reference_protein_ids)
    protein_ids_in_tax_level <- unique(c(monophyletic_clades_filter_selected_tax_level_sep_rows_curr$reference_protein_ids, unique_nonvertical_protein_ids))
    
    # Find most specific clade that includes all proteins in largest clade
    monophyletic_clades_filter_selected_tax_level_fraction_present <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% group_by(OG_id, label, clade_index) %>% summarize(fraction_parent_total_protein_ids_present_in_child = sum(protein_ids_in_tax_level %in% c(reference_protein_ids, unique(unlist(strsplit(nonvertical_protein_ids, split=","))))) / length(protein_ids_in_tax_level), fraction_parent_reference_protein_ids_present_in_child = sum(unique_reference_protein_ids %in% reference_protein_ids) / length(unique_reference_protein_ids), .groups = "drop")
    monophyletic_clades_filter_selected_tax_level_fraction_present <- monophyletic_clades_filter_selected_tax_level_fraction_present[order(monophyletic_clades_filter_selected_tax_level_fraction_present$fraction_parent_total_protein_ids_present_in_child, monophyletic_clades_filter_selected_tax_level_fraction_present$fraction_parent_reference_protein_ids_present_in_child, monophyletic_clades_filter_selected_tax_level_fraction_present$clade_index, decreasing=TRUE),]
    monophyletic_clades_filter_selected_tax_level_fraction_present <- monophyletic_clades_filter_selected_tax_level_fraction_present %>% filter(fraction_parent_total_protein_ids_present_in_child == 1)
    retain_clade_index <- monophyletic_clades_filter_selected_tax_level_fraction_present$clade_index[1]
    retain_clade_indexes <- c(retain_clade_indexes, retain_clade_index)
    
  }
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% filter(clade_index %in% retain_clade_indexes)
  monophyletic_clades_filter_selected_tax_level_sep_rows <- monophyletic_clades_filter_selected_tax_level_sep_rows[!duplicated(monophyletic_clades_filter_selected_tax_level_sep_rows$reference_protein_ids),]
  
  # If no rows remaining, skip
  if (nrow(monophyletic_clades_filter_selected_tax_level_sep_rows) == 0) {
    next
  }
  
  # Combine. Remove any nonvertical proteins not in selected tax level
  monophyletic_clades_filter_selected_tax_level_largest_clades <- monophyletic_clades_filter_selected_tax_level_sep_rows %>% group_by(OG_id, label, clade_index) %>% summarize(distance_to_root = mean(distance_to_root), reference_protein_ids = paste0(sort(unique(unlist(strsplit(reference_protein_ids, split=",")))), collapse=","), count = mean(count), species_overlap = mean(species_overlap), duplications_rec = mean(duplications_rec), nonvertical_protein_ids = paste0(sort(unique(unlist(strsplit(nonvertical_protein_ids, split=","))[gsub("_.*", "", unlist(strsplit(nonvertical_protein_ids, split=","))) %in% species_tree_subtree$tip.label])), collapse=","), n_species = length(unique(gsub("_.*", "", unlist(strsplit(reference_protein_ids, split=","))))), n_reference_proteins = length(unlist(strsplit(reference_protein_ids, split=","))), fraction_euk_species = mean(fraction_euk_species), mito_localization_prob = mean(mito_localization_prob), .groups = "drop")
  
  # Order by clade size from largest to smallest (by proteins, then species)
  monophyletic_clades_filter_selected_tax_level_largest_clades <- monophyletic_clades_filter_selected_tax_level_largest_clades[order(monophyletic_clades_filter_selected_tax_level_largest_clades$n_reference_proteins, monophyletic_clades_filter_selected_tax_level_largest_clades$n_species, decreasing=TRUE),]
  
  # Assign PhROG ids
  monophyletic_clades_filter_selected_tax_level_largest_clades$PhROG_id <- paste0(monophyletic_clades_filter_selected_tax_level_largest_clades$OG_id, "_", original_selected_tax_level, "_", sprintf("PhROG%07d", 1:nrow(monophyletic_clades_filter_selected_tax_level_largest_clades)))
  
  # Get vertical tax labels
  for (curr_rowindex in 1:nrow(monophyletic_clades_filter_selected_tax_level_largest_clades)) {
    curr_reference_protein_ids <- unlist(strsplit(monophyletic_clades_filter_selected_tax_level_largest_clades$reference_protein_ids[curr_rowindex], split=","))
    curr_reference_species_ids <- gsub("_.*", "", curr_reference_protein_ids)
    curr_reference_species_ids <- curr_reference_species_ids[curr_reference_species_ids %in% species_tree_subtree$tip.label]
    lca_index <- get_mrca_of_set(species_tree, curr_reference_species_ids)
    if (lca_index <= Ntip(species_tree)) {
      curr_species_tree_label <- species_tree$tip.label[lca_index]
    } else {
      curr_species_tree_label <- species_tree$node.label[lca_index - Ntip(species_tree)]
    }
    monophyletic_clades_filter_selected_tax_level_largest_clades$label[curr_rowindex] <- curr_species_tree_label
  }
  
  # Convert into long format.
  monophyletic_clades_filter_selected_tax_level_largest_clades_vertical_long <- monophyletic_clades_filter_selected_tax_level_largest_clades %>% separate_rows(reference_protein_ids, sep=",") %>% select(PhROG_id, reference_protein_ids, label, mito_localization_prob)
  monophyletic_clades_filter_selected_tax_level_largest_clades_nonvertical_long <- monophyletic_clades_filter_selected_tax_level_largest_clades %>% filter(nonvertical_protein_ids != "") %>% separate_rows(nonvertical_protein_ids, sep=",") %>% select(PhROG_id, nonvertical_protein_ids, label, mito_localization_prob)
  # Mark non-vertical proteins
  monophyletic_clades_filter_selected_tax_level_largest_clades_vertical_long$BOOL_NONVERTICAL <- FALSE
  monophyletic_clades_filter_selected_tax_level_largest_clades_nonvertical_long$BOOL_NONVERTICAL <- TRUE
  # Combine
  colnames(monophyletic_clades_filter_selected_tax_level_largest_clades_vertical_long)[2] <- c("protein_id")
  colnames(monophyletic_clades_filter_selected_tax_level_largest_clades_nonvertical_long)[2] <- c("protein_id")
  monophyletic_clades_filter_selected_tax_level_largest_clades_long <- rbind(monophyletic_clades_filter_selected_tax_level_largest_clades_vertical_long, monophyletic_clades_filter_selected_tax_level_largest_clades_nonvertical_long)
  monophyletic_clades_filter_selected_tax_level_largest_clades_long <- monophyletic_clades_filter_selected_tax_level_largest_clades_long[order(monophyletic_clades_filter_selected_tax_level_largest_clades_long$PhROG_id),]
  # Add primary OG label
  monophyletic_clades_filter_selected_tax_level_largest_clades_long$BOOL_primary_OG <- TRUE
  fusion_protein_primary_disjoint_OG_mapping_curr_OG <- fusion_protein_primary_disjoint_OG_mapping %>% filter(orthogroup_primary %in% monophyletic_clades_filter_selected_tax_level_largest_clades$OG_id)
  monophyletic_clades_filter_selected_tax_level_largest_clades_long$BOOL_primary_OG[which(monophyletic_clades_filter_selected_tax_level_largest_clades_long$protein_id %in% fusion_protein_primary_disjoint_OG_mapping$protein_id & !monophyletic_clades_filter_selected_tax_level_largest_clades_long$protein_id %in% fusion_protein_primary_disjoint_OG_mapping_curr_OG$protein_id)] <- FALSE
  monophyletic_clades_filter_selected_tax_level_largest_clades_primary_OG_summary <- monophyletic_clades_filter_selected_tax_level_largest_clades_long %>% group_by(PhROG_id) %>% filter(!BOOL_NONVERTICAL) %>% summarize(fraction_primary_OG_for_vertical_proteins = sum(BOOL_primary_OG) / n())
  monophyletic_clades_filter_selected_tax_level_largest_clades$fraction_primary_OG_for_vertical_proteins <- monophyletic_clades_filter_selected_tax_level_largest_clades_primary_OG_summary$fraction_primary_OG_for_vertical_proteins[match(monophyletic_clades_filter_selected_tax_level_largest_clades$PhROG_id, monophyletic_clades_filter_selected_tax_level_largest_clades_primary_OG_summary$PhROG_id)]
  
  # # Write out
  # write.table(monophyletic_clades_filter_selected_tax_level_largest_clades_long, paste0(curr_out_dir, "/", OG_id, "_", original_selected_tax_level, "_PhROGs_long.tsv"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
