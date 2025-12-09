# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(castor))

# Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
base_dir <- args[1]
split_transfer_trees_dir <- args[2]
out_dir <- args[3]
OG_id <- args[4]
method <- args[5]
BOOL_verbose <- FALSE

BOOL_split_transfers_small_clades <- TRUE
BOOL_filter_transfers_by_LCA_in_euk_supergroup <- TRUE
posttransfer_clade_species_threshold <- 10
posttransfer_clade_support_threshold <- 0.1
BOOL_largest_euk_monophyletic_clade <- FALSE
BOOL_exclude_prokaryotes <- FALSE
BOOL_soft_LCA <- FALSE
species_overlap_fraction_threshold <- 0.1

reference_species <- c("9606", "559292", "3702", "1257118", "5689", "185431", "5741", "36329", "32595", "508771")
incomplete_localization_data_species <- c("36329", "32595")

### Read in data
# Read in sample of reconciled trees
rec_trees_filename <- here("data/reconciled_trees", "/", OG_id, "/reconciliations/family_1.rec_uml") # for Alerax
trees <- read.tree(rec_trees_filename)

# Read in taxonomy data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
uniprot_proteomes_all_tax <- read.delim(here("data/taxonomy", "uniprot_all_tax_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), header=TRUE)
uniprot_proteomes_prok_tree_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain != "Eukaryota"]
uniprot_proteomes_euk_tree_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Eukaryota"]
uniprot_proteomes_bacteria_counts <- uniprot_proteomes_all_tax %>% filter(domain == "Bacteria") %>% group_by(tree_id) %>% summarize(n_species = n(), .groups="keep")
uniprot_proteomes_archaea_counts <- uniprot_proteomes_all_tax %>% filter(domain == "Archaea") %>% group_by(tree_id) %>% summarize(n_species = n(), .groups="keep")

# Read in species tree
species_tree <- here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample.contree")

# Calculate node depths
depths_df <- data.frame(label = c(species_tree$tip.label, species_tree$node.label), distance=get_all_distances_to_root(species_tree, as_edge_count=TRUE))

species_tree_labels <- c(species_tree$tip.label, species_tree$node.label)

# Generate long format taxonomy table based on species tree
# Start with tips
species_tree_tax <- data.frame(label = species_tree$tip.label, tree_id = species_tree$tip.label, depth = get_pairwise_distances(species_tree, species_tree$tip.label, rep(find_root(species_tree), Ntip(species_tree)), as_edge_counts = TRUE))
# Add in nodes
for (curr_node_label in species_tree$node.label) {
  curr_subtree <- get_subtree_at_node(species_tree, curr_node_label)$subtree
  depth <- get_pairwise_distances(species_tree, curr_node_label, find_root(species_tree), as_edge_counts = TRUE)
  species_tree_tax <- rbind(species_tree_tax, data.frame(label = curr_node_label, tree_id = curr_subtree$tip.label, depth = depth))
}
# Add index
species_tree_tax$index <- match(species_tree_tax$label, species_tree_labels)

# Get eukaryote subtree
euk_species_subtree <- drop.tip(species_tree, uniprot_proteomes_prok_tree_ids)
euk_species_subtree_labels <- c(euk_species_subtree$tip.label, euk_species_subtree$node.label)

# Get eukaryote proteins list
euk_protein_ids <- trees[[1]]$tip.label[gsub("_.*", "", trees[[1]]$tip.label) %in% uniprot_proteomes_euk_tree_ids]

# Determine which intergroup transfer events to consider that have LCA in a specified supergroup
permitted_split_transfer_donor_recipient_LCA_nodelabels <- unique(c("cellular_organisms", "Archaea", "Asgardgroup", "Eukaryota", "Opimoda", "Amorphea_CRuMs", "Amorphea", "Obazoa", "Diphoda", "Diaphorectickes", "CAM_Haptista", "CAM", "Archaeplastida", uniprot_proteomes_tax$superfamily[uniprot_proteomes_tax$domain == "Eukaryota"]))

# Read in the prokaryote mmseqs2 cluster data
prokaryote_mmseqs2_clusters_taxids <- read.table(here("data/downsample_prokaryotes", "prokaryote_mmseqs2_clusters_taxids.tsv"), sep="\t", quote="", header=TRUE)


## Read in mito-localization data
# Read in mito goldp gene list for binary labels: mito, non-mito
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2025.09.30.tsv"), sep="\t", header=TRUE)

# Read in DeepLoc results
deeploc_results <- read.table(file.path(supplemental_data_directory, "TableS7_retrained_DeepLoc_predictions.tsv"), header=TRUE)
colnames(deeploc_results) <- c("Protein_ID", "Mitochondrion")

# Read in organelle-encoded proteins
mtdna_proteins <- read.table(here("data/deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1
nonmito_organelle_proteins <- read.table(here("data/deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
missing_mtdna_taxids <- read.table(here("data/deeploc", "missing_mtdna_taxids.txt"))$V1
missing_nonmito_organelle_taxids <- read.table(here("data/deeploc", "missing_nonmito_organelle_taxids.txt"))$V1

# Ignore unedited mtDNA proteins and likely NUMTs in Arabidopsis
ath_proteins_curr <- grep("3702_", trees[[1]]$tip.label, fixed=TRUE, value=TRUE)
ath_uneditedmtdna_or_numt_proteins <- read.table(here("data/mito_orthogroups", "ath_unedited_numt_mtdna_proteins.txt"))$V1
ath_uneditedmtdna_or_numt_proteins <- ath_uneditedmtdna_or_numt_proteins[ath_uneditedmtdna_or_numt_proteins %in% ath_proteins_curr]
if (length(ath_uneditedmtdna_or_numt_proteins) > 0) {
  if (BOOL_verbose) {
    print(paste0("Ignoring Arabidopsis unedited mtDNA proteins: ", paste0(ath_uneditedmtdna_or_numt_proteins, collapse=",")))
  }
}

# Function to get immediate ancestor node
get_ancestor_node <- function(species_tree, index) {
  species_tree$edge[which(species_tree$edge[,2] == index),1]
}

### Identify and split recent horizontal gene transfers
if (BOOL_split_transfers_small_clades) {
  original_tree_index_list_filename <- paste0(split_transfer_trees_dir, "/split_transfer_trees/", OG_id, "_original_tree_index_list.txt")
  completed_trees_filename <- paste0(split_transfer_trees_dir, "/split_transfer_trees/", OG_id, "_split_transfer_trees.rec_uml")
  
  # Check to see if output already exists
  BOOL_split_transfers_small_clades_completed <- file.exists(completed_trees_filename)
  if (BOOL_split_transfers_small_clades_completed) {
    # If exists, read in
    original_tree_index_list <- read.table(original_tree_index_list_filename, sep="\t", header=FALSE)$V1
    completed_trees <- read.tree(completed_trees_filename)
  } else {
    # If doesn't exist, do work
    trees_to_process <- trees
    species_tree_raw <- species_tree
    
    if (method == "ALE") {
      species_tree_euks <- get_subtree_at_node(species_tree_raw, which(species_tree_raw$node.label == "1375"))$subtree # for ALE
    } else if (method == "AleRax") {
      species_tree_euks <- get_subtree_at_node(species_tree_raw, which(species_tree_raw$node.label == "Node34_Eukaryota"))$subtree # for Alerax
    }
    
    species_tree_proks <- drop.tip(species_tree_raw, species_tree_euks$tip.label)
    species_tree_euks_labels <- c(species_tree_euks$tip.label, species_tree_euks$node.label)
    species_tree_proks_labels <- c(species_tree_proks$tip.label, species_tree_proks$node.label)
    
    completed_trees <- c()
    original_tree_index_list <- 1:length(trees_to_process)
    original_tree_index_curr_index <- 1
    while(length(trees_to_process) > 0) {
      curr_tree <- trees_to_process[[1]]
      BOOL_updated_curr_tree <- FALSE
      
      if (BOOL_verbose) {
        print(length(trees_to_process))
        # print(curr_tree)
      }
      
      curr_tree_labels <- c(curr_tree$tip.label, curr_tree$node.label)
      
      transfer_node_labels_raw <- grep("->", curr_tree_labels, value=TRUE)
      transfer_node_index_raw <- grep("->", curr_tree_labels)
      node_labels_df <- data.frame(nodelabels=transfer_node_labels_raw, node_index=transfer_node_index_raw, transfer_node_full_labels = transfer_node_labels_raw)
      node_labels_df <- node_labels_df %>% separate_rows(nodelabels, sep="\\.")
      transfer_node_labels <- grep("T@", node_labels_df$nodelabels, fixed=TRUE, value=TRUE)
      transfer_node_indexes <- node_labels_df$node_index[grep("T@", node_labels_df$nodelabels, fixed=TRUE)]
      transfer_node_full_labels <- node_labels_df$transfer_node_full_labels[grep("T@", node_labels_df$nodelabels, fixed=TRUE)]
      
      if (Ntip(curr_tree) > 1 & length(transfer_node_labels) > 0) {
        for (curr_index in 1:length(transfer_node_labels)) {
          transfer_node_label <- transfer_node_labels[curr_index]
          transfer_node_full_label <- transfer_node_full_labels[curr_index]
          node_index <- transfer_node_indexes[curr_index]
          
          if (is.na(transfer_node_label)) {
            next
          }
          
          node_labels_df <- data.frame(nodelabels=transfer_node_label)
          node_labels_df <- node_labels_df %>% separate_rows(nodelabels, sep="\\.")
          transfer_node_id_raw <- grep("T@", node_labels_df$nodelabels, fixed=TRUE, value=TRUE)
          transfer_event_id <- gsub(".*(T@[a-zA-Z0-9_]+->[a-zA-Z0-9_]+).*", "\\1", transfer_node_id_raw)
          transfer_event_donor <- gsub("T@", "", gsub("->.*", "", transfer_event_id))
          transfer_event_recipient <- unique(gsub(".*->", "", transfer_event_id))
          
          if (BOOL_verbose) {
            print(transfer_node_label)
            print(transfer_node_full_label)
            print(transfer_node_id_raw)
            print(transfer_event_id)
            print(transfer_event_recipient)
          }
          
          if (length(transfer_event_recipient) > 1) {
            print(paste0("Warning: Length of transfer_event_recipient is more than 1 for ", transfer_event_recipient, " for ", transfer_node_label, " with recipient ", transfer_event_recipient))
            next
          }
          
          # Get tip/node index of daughter with the recipient clade
          daughter_indexes <- curr_tree$edge[which(curr_tree$edge[,1] == node_index),2]
          daughter_node_labels <- c()
          subtrees <- c()
          for (daughter_index in daughter_indexes) {
            if (daughter_index > Ntip(curr_tree)) {
              # If it's a node
              daughter_node_label <- curr_tree$node.label[daughter_index - Ntip(curr_tree)]
              # daughter_node_label <- gsub(".*\\.", "", daughter_node_label) # Get just the oldest event at each node
              daughter_node_labels <- c(daughter_node_labels, daughter_node_label)
              if (length(subtrees) == 0) {
                subtrees <- c(get_subtree_at_node(curr_tree, daughter_index - Ntip(curr_tree))$subtree)
              } else {
                subtrees <- c(subtrees, get_subtree_at_node(curr_tree, daughter_index - Ntip(curr_tree))$subtree)
              }
            } else {
              # If it's a tip, extract the taxid
              # daughter_node_label <- curr_tree$tip.label[daughter_index]
              daughter_node_label <- gsub("_.*", "", curr_tree$tip.label[daughter_index])
              daughter_node_labels <- c(daughter_node_labels, daughter_node_label)
              if (length(subtrees) == 0) {
                subtrees <- c(keep.tip(curr_tree, curr_tree$tip.label[daughter_index]))
              } else {
                subtrees <- c(subtrees, keep.tip(curr_tree, curr_tree$tip.label[daughter_index]))
              }
            }
          }
          
          # If either daughter node labels has a transfer, skip
          if (any(grep("T@", daughter_node_labels, fixed=TRUE))) {
            next
          }
          # Assign subtrees
          subtree_1 <- subtrees[[1]]
          subtree_2 <- subtrees[[2]]
          
          # If sister subtrees are in different supergroups, have sufficient phylogenetic discordance to infer a possible transfer
          mrca_index <- get_mrca_of_set(species_tree, c(gsub("_.*", "", subtree_1$tip.label), gsub("_.*", "", subtree_2$tip.label)))
          mrca_node_label <- species_tree_labels[mrca_index]
          mrca_tax_label <- gsub("Node[0-9]+_", "", mrca_node_label)

          if (!BOOL_filter_transfers_by_LCA_in_euk_supergroup | mrca_tax_label %in% permitted_split_transfer_donor_recipient_LCA_nodelabels) {
            ## Determine direction of transfer
            # Get cousin clade
            parent_index <- curr_tree$edge[curr_tree$edge[,2] == node_index,1]
            if (length(parent_index) == 0) {
              next
            }
            pretransfer_subtree <- get_subtree_at_node(curr_tree, parent_index - Ntip(curr_tree))$subtree
            
            # Drop either subtree
            cousin_drop_subtree1 <- drop.tip(pretransfer_subtree, c(subtree_1$tip.label))
            cousin_drop_subtree2 <- drop.tip(pretransfer_subtree, c(subtree_2$tip.label))
            cousin_drop_subtree1_mrca_index <- get_mrca_of_set(species_tree, c(gsub("_.*", "", cousin_drop_subtree1$tip.label)))
            cousin_drop_subtree1_mrca_node_label <- species_tree_labels[cousin_drop_subtree1_mrca_index]
            cousin_drop_subtree2_mrca_index <- get_mrca_of_set(species_tree, c(gsub("_.*", "", cousin_drop_subtree2$tip.label)))
            cousin_drop_subtree2_mrca_node_label <- species_tree_labels[cousin_drop_subtree2_mrca_index]
            
            dist_to_root_cousin_drop_subtree1 <- get_pairwise_distances(species_tree, find_root(species_tree), cousin_drop_subtree1_mrca_node_label, as_edge_counts = TRUE)
            dist_to_root_cousin_drop_subtree2 <- get_pairwise_distances(species_tree, find_root(species_tree), cousin_drop_subtree2_mrca_node_label, as_edge_counts = TRUE)
            
            if (BOOL_verbose) {
              print(transfer_event_donor)
              print(transfer_event_recipient)
              print(cousin_drop_subtree1_mrca_node_label)
              print(cousin_drop_subtree2_mrca_node_label)
            }
            
            # Determine direction of transfer by vertical similarity: prefer more recent monophyletic clade
            if (dist_to_root_cousin_drop_subtree1 > dist_to_root_cousin_drop_subtree2) {
              # Replace Alerax recipient with the vertical donor
              subtree_to_split <- subtree_1
              not_subtree_to_split <- subtree_2
              vertical_donor <- cousin_drop_subtree2_mrca_node_label
            } else if (dist_to_root_cousin_drop_subtree1 < dist_to_root_cousin_drop_subtree2) {
              # Alerax recipient agrees with vertical donor
              subtree_to_split <- subtree_2
              not_subtree_to_split <- subtree_1
              vertical_donor <- cousin_drop_subtree1_mrca_node_label
            } else {
              if (BOOL_verbose) {
                print("Ambiguous donor/recipient")
              }
              next
            }
            
            # Only split clades with few species. Take mmseqs2 downsampled prokaryotes into account.
            curr_tree_ids <- gsub("_.*", "", subtree_to_split$tip.label)
            curr_tree_ids_euk <- curr_tree_ids[which(curr_tree_ids %in% uniprot_proteomes_euk_tree_ids)]
            proteins_remove_first_underscore <- sub("^[^_]*_", "", subtree_to_split$tip.label)
            prok_proteins_remove_first_underscore <- proteins_remove_first_underscore[which(curr_tree_ids %in% uniprot_proteomes_prok_tree_ids)]
            if (length(prok_proteins_remove_first_underscore) > 0) {
              prokaryote_mmseqs2_clusters_taxids_in_tree <- prokaryote_mmseqs2_clusters_taxids[which(prokaryote_mmseqs2_clusters_taxids$rep_seq %in% prok_proteins_remove_first_underscore),]
              prok_species_in_clade <- unlist(lapply(prok_proteins_remove_first_underscore, function(lst) {
                unique(unlist(strsplit(prokaryote_mmseqs2_clusters_taxids_in_tree$unique_member_seq_taxids[which(prokaryote_mmseqs2_clusters_taxids_in_tree$rep_seq %in% lst)], split=",")))
              }))
              n_unique_species_subtree_to_split <- length(unique(unlist(prok_species_in_clade))) + length(unique(curr_tree_ids_euk))
            } else {
              n_unique_species_subtree_to_split <- length(unique(curr_tree_ids_euk))
            }
            
            if (n_unique_species_subtree_to_split > posttransfer_clade_species_threshold) {
              next
            }
            
            if (BOOL_verbose) {
              vertical_recipient <- species_tree_labels[get_mrca_of_set(species_tree, c(gsub("_.*", "", subtree_to_split$tip.label)))]
              print(paste0("Vertically inferred transfer: ", vertical_donor, "->", vertical_recipient))
              print(not_subtree_to_split)
              print(subtree_to_split)
            }
            
            # Update trees
            trees_to_process[[1]] <- drop.tip(curr_tree, subtree_to_split$tip.label)
            trees_to_process <- c(trees_to_process, subtree_to_split)
            original_tree_index_list <- c(original_tree_index_list, original_tree_index_list[original_tree_index_curr_index])
            BOOL_updated_curr_tree <- TRUE
            break
          }
        }
      }
      
      # If current tree hasn't been updated, update the completed and remaining tree set
      if (!BOOL_updated_curr_tree) {
        # Add to completed tree set
        if (length(completed_trees) == 0) {
          completed_trees <- c(curr_tree)
        } else {
          completed_trees <- c(completed_trees, curr_tree)
        }
        if (BOOL_verbose) {
          print(paste0("Processing completed for ", length(completed_trees)))
        }
        original_tree_index_curr_index <- original_tree_index_curr_index + 1
        
        # Update set of trees to be processed
        if (length(trees_to_process) > 1) {
          trees_to_process <- trees_to_process[-1]
        } else {
          trees_to_process <- NULL # no more trees
        }
      }
    }
    
    ## Write out the split transfer trees
    # write.table(original_tree_index_list, original_tree_index_list_filename, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    # write.tree(completed_trees, completed_trees_filename)
  }
}

### Get full set of posterior clades. For each clade, identify duplications and infer crude ancestral mito localization
vertical_result_agg <- c()

for (i in 1:length(trees)) {
  if (BOOL_verbose) {
    print(paste0(i, " of ", length(trees)))
  }
  
  curr_tree <- trees[[i]]
  
  curr_tree_mito_localization_states <- rep(NA, Ntip(curr_tree))
  
  protein_accessions <- curr_tree$tip.label
  
  # Remove transcript id for ACANB to match the goldp protein ids
  if (any(grepl("^1257118_", protein_accessions))) {
    protein_accessions_reformat <- protein_accessions
    protein_accessions_reformat[grep("^1257118_", protein_accessions_reformat)] <- gsub("_t.*", "", protein_accessions_reformat[grep("^1257118_", protein_accessions_reformat)])
  } else {
    protein_accessions_reformat <- protein_accessions
  }
  
  # Mark states for reference species with complete localization data: binary mito or not mito
  curr_tree_mito_localization_states[which(gsub("_.*", "", curr_tree$tip.label) %in% reference_species[!reference_species %in% incomplete_localization_data_species])] <- 1 # nonmito
  curr_tree_mito_localization_states[which(protein_accessions_reformat %in% gold_gene_accession_OG_id_df$gene_accession)] <- 2 # mito
  
  # Ignore unedited mtDNA proteins and likely NUMTs in Arabidopsis for localization
  if (length(ath_uneditedmtdna_or_numt_proteins) > 0) {
    curr_tree_mito_localization_states[which(protein_accessions_reformat %in% ath_uneditedmtdna_or_numt_proteins)] <- NA # ignore
  }
  
  ## Use prior
  default_prior <- NA # no prior information for unannotated species
  tip_priors_mat <- matrix(data = default_prior, nrow = length(curr_tree_mito_localization_states), ncol = 2)
  
  # Mark deeploc predictions as priors
  tip_priors_mat[which(protein_accessions %in% deeploc_results$Protein_ID),2] <- deeploc_results$Mitochondrion[match(protein_accessions[which(protein_accessions %in% deeploc_results$Protein_ID)], deeploc_results$Protein_ID)]
  
  # Mark organelle-encoded proteins. Omit deeploc predictions for species without mtDNA annotation
  n_mtdna_proteins <- sum(protein_accessions %in% mtdna_proteins)
  n_nonmito_organelle_proteins <- sum(protein_accessions %in% nonmito_organelle_proteins)
  if (n_mtdna_proteins > 0) {
    # Ignore deeploc predictions for species missing mtDNA annotation (excepting species where no mtDNA exists)
    omit_deeploc_index <- which(protein_accessions %in% deeploc_results$Protein_ID & gsub("_.*", "", protein_accessions) %in% missing_mtdna_taxids)
    tip_priors_mat[omit_deeploc_index,2] <- NA
  }
  if (n_nonmito_organelle_proteins > 0) {
    # Ignore deeploc predictions for species missing mtDNA annotation (excepting species where no mtDNA exists)
    omit_deeploc_index <- which(protein_accessions %in% deeploc_results$Protein_ID & gsub("_.*", "", protein_accessions) %in% missing_nonmito_organelle_taxids)
    tip_priors_mat[omit_deeploc_index,2] <- NA
  }
  tip_priors_mat[which(protein_accessions %in% nonmito_organelle_proteins),2] <- 0
  tip_priors_mat[which(protein_accessions %in% mtdna_proteins),2] <- 1
  
  # Mark known mito proteins as certain. Omit nonmito proteins since some may be truly mito but not captured by proteomics
  tip_priors_mat[which(curr_tree_mito_localization_states == 2),2] <- 1
  
  # Mark nonmito state as the complement of mito state
  tip_priors_mat[,1] <- 1 - tip_priors_mat[,2]
  
  if (sum(!is.na(tip_priors_mat[,2])) > 1) {
    # If there are multiple scores available
    result <- hsp_mk_model(tree=curr_tree, tip_states=NULL, tip_priors=tip_priors_mat, Nstates=2, rate_model = "ER", Ntrials=5)
    
    mito_likelihood <- result$likelihoods[,2] # take the mito state likelihood
    
  } else if (sum(!is.na(tip_priors_mat[,2])) == 1) {
    # If there is only 1 protein with a score, just use its score to inpute the rest
    mito_likelihood <- rep(tip_priors_mat[!is.na(tip_priors_mat[,2]),2], Ntip(curr_tree) + Nnode(curr_tree))
    
  } else {
    # If there are no scores available, assume nonmito and set to 0
    print(paste0("No annotated species proteins present in ", OG_id))
    mito_likelihood <- rep(0, Ntip(curr_tree) + Nnode(curr_tree))
  }
  
  # Get species tree labels, descendant protein ids, and duplications by species overlap.
  species_tree_label_vertical <- c()
  reference_protein_ids_vertical <- c()
  reference_protein_ids_nonvertical <- c()
  species_overlap_vertical <- c()
  species_overlap_vertical_binary <- c()
  species_overlap_taxids <- c()
  duplications_rec <- c()
  for (node_index in 1:Nnode(curr_tree)) {
    curr_subtree <- get_subtree_at_node(curr_tree, node_index)$subtree
    
    curr_tree_ids <- unique(gsub("_.*", "", curr_subtree$tip.label))
    lca_index <- get_mrca_of_set(species_tree, curr_tree_ids)
    
    if (lca_index <= Ntip(species_tree)) {
      curr_species_tree_label <- species_tree$tip.label[lca_index]
    } else {
      curr_species_tree_label <- species_tree$node.label[lca_index - Ntip(species_tree)]
    }
    
    species_tree_label_vertical <- c(species_tree_label_vertical, curr_species_tree_label)

    # Compute species overlap
    if (Ntip(curr_subtree) > 1) {
      descendant_indexes <- curr_subtree$edge[curr_subtree$edge[,1] == find_root(curr_subtree),2]
      
      descendant_tip_labels <- c()
      for (descendant_index in descendant_indexes) {
        if (descendant_index <= Ntip(curr_subtree)) {
          curr_descendant_tip_labels <- curr_subtree$tip.label[descendant_index]
        } else {
          curr_descendant_tip_labels <- get_subtree_at_node(curr_subtree, descendant_index - Ntip(curr_subtree))$subtree$tip.label
        }
        descendant_tip_labels <- c(descendant_tip_labels, list(curr_descendant_tip_labels))
      }
      
      descendant_species <- lapply(descendant_tip_labels, function(lst) unique(gsub("_.*", "", unlist(strsplit(lst, split=",")))))
      descendant_species <- lapply(descendant_species, function(lst) lst[lst %in% uniprot_proteomes_euk_tree_ids]) # filter for euks
      
      if (length(descendant_species[[1]]) > 0 & length(descendant_species[[2]]) > 0) {
        # Use fraction species overlap
        species_overlap_fraction <- sum(duplicated(unlist(descendant_species))) / length(unique(unlist(descendant_species)))
        species_overlap_vertical <- c(species_overlap_vertical, species_overlap_fraction)
        curr_taxids_duplicated <- sort(unique(unlist(descendant_species)[duplicated(unlist(descendant_species))]))
        species_overlap_taxids <- c(species_overlap_taxids, paste0(curr_taxids_duplicated, collapse=","))
        
        # Use binarized fraction species overlap
        species_overlap_fraction <- sum(duplicated(unlist(descendant_species))) / max(1, length(unique(unlist(descendant_species))))
        if (species_overlap_fraction > species_overlap_fraction_threshold) {
          species_overlap_vertical_binary <- c(species_overlap_vertical_binary, 1)
        } else {
          species_overlap_vertical_binary <- c(species_overlap_vertical_binary, 0)
        }
        
      } else {
        species_overlap_vertical <- c(species_overlap_vertical, 0)
        species_overlap_vertical_binary <- c(species_overlap_vertical_binary, 0)
        species_overlap_taxids <- c(species_overlap_taxids, "")
      }
    } else {
      species_overlap_vertical <- c(species_overlap_vertical, 0)
      species_overlap_vertical_binary <- c(species_overlap_vertical_binary, 0)
      species_overlap_taxids <- c(species_overlap_taxids, "")
    }
    
    # Get number of duplications at the current root node
    duplications_rec <- c(duplications_rec, length(grep("D@", curr_subtree$node.label[1], fixed=TRUE)))
    
    # Store reference species protein ids. Use the raw nonsplit transfer subtrees, otherwise large clades become highly fragmented.
    reference_protein_ids_vertical <- c(reference_protein_ids_vertical, paste0(curr_subtree$tip.label, collapse=","))
    
    # Get non-vertical protein ids
    if (BOOL_split_transfers_small_clades) {
      # Get LCA from the split transfer trees
      curr_subtree_node <- curr_tree$node.label[node_index]
      original_tree_index_list_curr <- which(original_tree_index_list == i)
      curr_subtree_split_transfers_tiplabels <- c()
      fraction_proteins_covered <- c()
      count <- 0
      for (j in 1:length(original_tree_index_list_curr)) {
        original_tree_index_curr <- original_tree_index_list_curr[j]
        
        if (curr_subtree_node %in% completed_trees[[original_tree_index_curr]]$node.label) {
          curr_subtree_split_transfers <- get_subtree_at_node(completed_trees[[original_tree_index_curr]], curr_subtree_node)$subtree
          if (all(curr_subtree_split_transfers$tip.label %in% curr_subtree$tip.label)) {
            curr_subtree_split_transfers_tiplabels <- c(curr_subtree_split_transfers_tiplabels, curr_subtree_split_transfers$tip.label)
            count <- count + 1
          }
        }
        fraction_proteins_covered <- c(fraction_proteins_covered, sum(curr_subtree$tip.label %in% completed_trees[[original_tree_index_curr]]$tip.label) / length(curr_subtree$tip.label))
      }
      # If have no matching nodelabels, identify the donor tree as the oldest tree that contains any of the current subtree proteins.
      if (count == 0) {
        parent_tree_index <- original_tree_index_list_curr[which(fraction_proteins_covered > 0)[1]]
        curr_subtree_split_transfers_tiplabels <- completed_trees[[parent_tree_index]]$tip.label
      }
      
      reference_protein_ids_nonvertical <- c(reference_protein_ids_nonvertical, paste0(sort(curr_subtree$tip.label[!curr_subtree$tip.label %in% curr_subtree_split_transfers_tiplabels]), collapse=","))
      
    } else {
      reference_protein_ids_nonvertical <- c(reference_protein_ids_nonvertical, "")
    }
  }
  
  reference_protein_accessions <- protein_accessions

  curr_result_df <- data.frame(OG_id=OG_id, label_vertical=c(gsub("_.*", "", protein_accessions), species_tree_label_vertical), reference_protein_ids=c(reference_protein_accessions, reference_protein_ids_vertical), reference_protein_ids_nonvertical = c(rep("", length(protein_accessions)), reference_protein_ids_nonvertical), species_overlap = c(rep(0, length(protein_accessions)), species_overlap_vertical), species_overlap_binary = c(rep(0, length(protein_accessions)), species_overlap_vertical_binary), species_overlap_taxids = c(rep("", length(protein_accessions)), species_overlap_taxids), duplications_rec = c(rep(0, length(protein_accessions)), duplications_rec), mito_localization_prob = mito_likelihood)
  curr_result_df$distance_to_root <- depths_df$distance[match(curr_result_df$label_vertical, depths_df$label)]
  
  # Split reconciliations
  curr_result_df_sep_rows <- curr_result_df
  curr_result_df_sep_rows$clade_index <- 1:nrow(curr_result_df_sep_rows)
  
  # Keep largest euk monophyletic clades
  if (BOOL_largest_euk_monophyletic_clade) {
    curr_result_df_sep_rows <- curr_result_df_sep_rows %>% filter(label_vertical %in% c(euk_species_subtree_labels, euk_protein_ids))
    curr_result_df_sep_rows <- curr_result_df_sep_rows %>% rowwise() %>% mutate(n_proteins = length(unique(unlist(strsplit(reference_protein_ids, split=",")))))
    curr_result_df_sep_rows <- curr_result_df_sep_rows[order(c(curr_result_df_sep_rows$n_proteins, 1/curr_result_df_sep_rows$distance_to_root), decreasing=TRUE),]
    curr_result_df_sep_rows <- curr_result_df_sep_rows %>% separate_rows(reference_protein_ids, sep=",")
    curr_result_df_sep_rows <- curr_result_df_sep_rows[!duplicated(curr_result_df_sep_rows$reference_protein_ids),]
  }
  
  vertical_result <- curr_result_df_sep_rows %>% group_by(OG_id, label_vertical, clade_index) %>% summarize(distance_to_root = mean(distance_to_root), reference_protein_ids = paste0(sort(unique(unlist(strsplit(reference_protein_ids, split=",")))), collapse=","), reference_protein_ids_nonvertical = paste0(sort(unique(unlist(strsplit(reference_protein_ids_nonvertical, split=",")))), collapse=","), species_overlap = mean(species_overlap), species_overlap_binary = mean(species_overlap_binary), species_overlap_taxids = paste0(sort(unique(unlist(strsplit(species_overlap_taxids, split=",")))), collapse=","), duplications_rec = mean(duplications_rec), mito_localization_prob = mean(mito_localization_prob), .groups = "keep")
  
  # Remove duplicates that can occur when dropping non-vertical proteins
  vertical_result <- vertical_result[!duplicated(vertical_result[,c("label_vertical", "reference_protein_ids")]),]
  
  # When using reference species protein ids to define groups
  vertical_result_agg <- rbind(vertical_result_agg, vertical_result)
}

### Remove well-supported non-vertical proteins and update ancestral species labels

## Ignore well-supported non-vertical proteins
if (BOOL_split_transfers_small_clades) {
  if (BOOL_verbose) {
    print("Finding well-supported non-vertical proteins")
  }
  # Find well-supported non-vertical proteins, at the per protein level
  vertical_result_summary_count_per_clade <- vertical_result_agg %>% group_by(OG_id, label_vertical, distance_to_root, reference_protein_ids) %>% summarize(count_per_clade = n(), .groups = "drop")
  vertical_result_summary_count_per_nonvertical <- vertical_result_agg %>% separate_rows(reference_protein_ids_nonvertical, sep=",") %>% filter(reference_protein_ids_nonvertical != "") %>% group_by(OG_id, label_vertical, distance_to_root, reference_protein_ids, reference_protein_ids_nonvertical) %>% summarize(count_per_nonvertical_protein = n(), .groups = "drop")
  vertical_result_summary_count_per_nonvertical$count_per_clade <- vertical_result_summary_count_per_clade$count_per_clade[match(vertical_result_summary_count_per_nonvertical$reference_protein_ids, vertical_result_summary_count_per_clade$reference_protein_ids)]
  
  vertical_result_summary_count_per_nonvertical_majority <- vertical_result_summary_count_per_nonvertical %>% mutate(fraction_nonvertical_protein_per_clade = count_per_nonvertical_protein / count_per_clade)

  # Require that a protein is non-vertical in a certain fraction of the clade in order to assign as non-vertical
  vertical_result_summary_count_per_nonvertical_majority <- vertical_result_summary_count_per_nonvertical_majority %>% filter(fraction_nonvertical_protein_per_clade >= posttransfer_clade_support_threshold)
  vertical_result_summary_count_per_nonvertical_majority_summary <- vertical_result_summary_count_per_nonvertical_majority %>% group_by(OG_id, label_vertical, distance_to_root, reference_protein_ids) %>% summarize(reference_protein_ids_nonvertical_majority = paste0(sort(unique(reference_protein_ids_nonvertical)), collapse=","), .groups = "drop")
  vertical_result_agg$reference_protein_ids_nonvertical <- vertical_result_summary_count_per_nonvertical_majority_summary$reference_protein_ids_nonvertical_majority[match(vertical_result_agg$reference_protein_ids, vertical_result_summary_count_per_nonvertical_majority_summary$reference_protein_ids)]
  vertical_result_agg$reference_protein_ids_nonvertical[is.na(vertical_result_agg$reference_protein_ids_nonvertical)] <- ""
} else {
  vertical_result_agg$reference_protein_ids_nonvertical <- ""
}

vertical_result_summary <- vertical_result_agg %>% group_by(OG_id, label_vertical, distance_to_root, reference_protein_ids, reference_protein_ids_nonvertical) %>% summarize(count = n(), species_overlap = mean(species_overlap), species_overlap_binary = mean(species_overlap_binary), species_overlap_taxids = paste0(sort(unique(unlist(strsplit(species_overlap_taxids, split=",")))), collapse=","), duplications_rec = mean(duplications_rec), mito_localization_fraction = mean(mito_localization_prob), .groups = "keep")
vertical_result_summary <- vertical_result_summary[!is.na(vertical_result_summary$label_vertical),]
protein_mito_likelihood_summary <- protein_mito_likelihood_agg %>% group_by(OG_id, protein_id) %>% summarize(mito_localization_fraction = mean(mito_localization_prob), .groups = "keep")


## Updating ancestral species labels and removing non-vertical proteins
if (BOOL_split_transfers_small_clades) {
  # Update ancestral species labels
  for (i in 1:nrow(vertical_result_summary)) {
    if (vertical_result_summary$reference_protein_ids_nonvertical[i] == "") {
      next
    }
    
    curr_reference_protein_ids <- unlist(strsplit(vertical_result_summary$reference_protein_ids[i], split=","))
    curr_reference_protein_ids_nonvertical <- unlist(strsplit(vertical_result_summary$reference_protein_ids_nonvertical[i], split=","))
    curr_reference_protein_ids_exclude_nonvertical <- curr_reference_protein_ids[!curr_reference_protein_ids %in% curr_reference_protein_ids_nonvertical]
    curr_reference_species_ids_exclude_nonvertical <- gsub("_.*", "", curr_reference_protein_ids_exclude_nonvertical)
    lca_index <- get_mrca_of_set(species_tree, curr_reference_species_ids_exclude_nonvertical)
    if (lca_index <= Ntip(species_tree)) {
      curr_species_tree_label <- species_tree$tip.label[lca_index]
    } else {
      curr_species_tree_label <- species_tree$node.label[lca_index - Ntip(species_tree)]
    }
    
    vertical_result_summary$label_vertical[i] <- curr_species_tree_label
    vertical_result_summary$reference_protein_ids[i] <- paste0(sort(curr_reference_protein_ids_exclude_nonvertical), collapse=",")
  }
}

vertical_result_summary <- vertical_result_summary %>% rowwise() %>% mutate(n_species = length(unique(gsub("_.*", "", unlist(strsplit(reference_protein_ids, split=","))))), n_reference_proteins = length(unique(unlist(strsplit(reference_protein_ids, split=",")))))
vertical_result_summary <- vertical_result_summary[order(vertical_result_summary$n_reference_proteins, 1/vertical_result_summary$distance_to_root, vertical_result_summary$count, decreasing=TRUE),]

# Write out
# write.table(vertical_result_summary, paste0(out_dir, "/", OG_id, "_euk_monophyletic_clades.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



