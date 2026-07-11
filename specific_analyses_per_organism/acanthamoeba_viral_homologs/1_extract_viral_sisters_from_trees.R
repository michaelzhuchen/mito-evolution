
# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(phytools))
suppressMessages(library(castor))

treefile_list <- "ACANB_virus_alignments_and_trees_highest_LL_trees.txt"
treefile_dir <- "ACANB_virus_alignments_and_trees_highest_LL_trees"

# Set results directory
result_directory <- here("specific_analyses_per_organism", "acanthamoeba_viral_homologs")

# Read in uniprot tax data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
euk_taxids <- uniprot_proteomes_tax$TaxId[which(uniprot_proteomes_tax$domain == "Eukaryota")]

# Virus taxids
virus_taxids <- c("212035", "554168", "2080449")

## Select donor/recipient taxa
donor_name <- "virus"
donor_taxids <- virus_taxids
recipient_taxids <- "1257118"

clade_purity_threshold <- 0.5

treefiles <- read.table(treefile_list, sep="\t", header=FALSE)$V1
treefiles <- sort(treefiles)

hgt_results_agg <- c()

for (treefile in treefiles) {
  print(treefile)
  
  tree <- read.tree(paste0(treefile_dir, "/", treefile))
  tree_labels <- c(tree$tip.label, tree$node.label)
  
  recipient_indexes <- which(gsub("_.*", "", tree$tip.label) %in% recipient_taxids)
  
  # Get immediate ancestor node
  get_ancestor_node <- function(tree, index) {
    tree$edge[which(tree$edge[,2] == index),1]
  }
  
  for (recipient_index in recipient_indexes) {
    recipient_protein <- tree$tip.label[recipient_index]
    
    farthest_tips <- find_farthest_tips(tree, only_descending_tips = FALSE, target_tips = NULL, as_edge_counts = TRUE, check_input = TRUE)
    curr_farthest_tip <- farthest_tips$farthest_tip_per_tip[recipient_index]
    
    curr_tree <- root(tree, outgroup = curr_farthest_tip, resolve.root=TRUE)
    
    # Climb up the tree until finding a clade containing < threshold recipient proteins, to designate as the recipient clade
    sister_clade_edge_index <- get_ancestor_node(curr_tree, recipient_index)
    sister_subtree <- get_subtree_at_node(curr_tree, sister_clade_edge_index - Ntip(curr_tree))$subtree
    other_sister_tiplabels <- sister_subtree$tip.label[!gsub("_.*", "", sister_subtree$tip.label) %in% recipient_taxids]
    recipient_sister_tiplabels <- sister_subtree$tip.label[gsub("_.*", "", sister_subtree$tip.label) %in% recipient_taxids]
    sister_node_support <- sister_subtree$node.label[1]
    while (length(recipient_sister_tiplabels) / length(sister_subtree$tip.label) >= clade_purity_threshold) {
      sister_clade_edge_index <- get_ancestor_node(curr_tree, sister_clade_edge_index)
      if (sister_clade_edge_index == find_root(curr_tree)) {
        sister_subtree <- curr_tree
        other_sister_tiplabels <- sister_subtree$tip.label[!gsub("_.*", "", sister_subtree$tip.label) %in% recipient_taxids]
        recipient_sister_tiplabels <- sister_subtree$tip.label[gsub("_.*", "", sister_subtree$tip.label) %in% recipient_taxids]
        sister_node_support <- NA
        break
      } else {
        sister_subtree <- get_subtree_at_node(curr_tree, sister_clade_edge_index - Ntip(curr_tree))$subtree
        other_sister_tiplabels <- sister_subtree$tip.label[!gsub("_.*", "", sister_subtree$tip.label) %in% recipient_taxids]
        recipient_sister_tiplabels <- sister_subtree$tip.label[gsub("_.*", "", sister_subtree$tip.label) %in% recipient_taxids]
        sister_node_support <- sister_subtree$node.label[1]
      }
    }
    
    cousin_clade_edge_index <- curr_tree$edge[curr_tree$edge[,2] == sister_clade_edge_index,1]
    if (sister_clade_edge_index != find_root(curr_tree)) {
      cousin_subtree <- get_subtree_at_node(curr_tree, cousin_clade_edge_index - Ntip(curr_tree))$subtree
      other_cousin_tiplabels <- cousin_subtree$tip.label[!gsub("_.*", "", cousin_subtree$tip.label) %in% recipient_taxids]
      other_cousin_tiplabels <- other_cousin_tiplabels[!other_cousin_tiplabels %in% sister_subtree$tip.label]
      cousin_node_support <- cousin_subtree$node.label[1]
    } else {
      other_cousin_tiplabels <- NA
      cousin_node_support <- NA
    }

    # Filter for donor proteins
    other_sister_tiplabels_donors <- other_sister_tiplabels[gsub("_.*", "", other_sister_tiplabels) %in% donor_taxids]
    other_cousin_tiplabels_donors <- other_cousin_tiplabels[gsub("_.*", "", other_cousin_tiplabels) %in% donor_taxids]

    hgt_results <- data.frame(Orthogroup = gsub("\\..*", "", treefile), recipient_protein, selected_donor_sister_tiplabels = paste0(sort(unique(other_sister_tiplabels_donors)), collapse=","), selected_donor_cousin_tiplabels = paste0(sort(unique(other_cousin_tiplabels_donors)), collapse=","), all_sister_tiplabels = paste0(sort(unique(other_sister_tiplabels)), collapse=","), all_cousin_tiplabels = paste0(sort(unique(other_cousin_tiplabels)), collapse=","), sister_node_support, cousin_node_support)
    hgt_results_agg <- rbind(hgt_results_agg, hgt_results)
  }
}

## Write out
# write.table(hgt_results_agg, paste0(result_directory, "/", "HGT_results_", donor_name, "_to_", paste0(recipient_taxids, collapse="."), ".tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote = TRUE)

