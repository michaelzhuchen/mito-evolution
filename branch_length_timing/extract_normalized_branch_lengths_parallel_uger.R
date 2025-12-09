# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(castor))
suppressMessages(library(phytools))

# Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
base_dir <- args[1]
phrogs_dir <- args[2]
alerax_dir <- args[3]
out_dir <- args[4]
bool_exclude_nonvertical <- TRUE
bool_exclude_small_clades <- TRUE
bool_skip_completed <- TRUE
OG_list <- args[8]
selected_tax_level <- args[9]
bool_verbose <- FALSE

## Read in data
# Read in species tree for euk203spp_prokgroups
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

n_minimum_species_leca <- 4


## Read in PhROGs
# Read in Eukaryota parent PhROGs
phrogs_long_eukaryota_parent <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long_eukaryota_parent) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_parent <- phrogs_long_eukaryota_parent %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))
# Read in LECA PhROGs
phrogs_long_eukaryota <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long_eukaryota) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota <- phrogs_long_eukaryota %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))

# Read in Lta mapping for Lta2019 to Ltaref
lta_mapping <- read.delim(here("data/orthogroups/idmapping", "map.ltaref.lta.exact.txt"), header=FALSE)
colnames(lta_mapping) <- c("ltaref_id", "lta2019_id")
lta_mapping$ltaref_id <- paste0("5689_", lta_mapping$ltaref_id)
lta_mapping$lta2019_id[lta_mapping$lta2019_id == ""] <- NA
lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)] <- paste0("5689_", lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
# Edge case: Drop one-to-many lta2019->ltaref mappings to mismatched proteins between tree and PhROGs
lta2019_onetomany_ids <- lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)][which(duplicated(lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)]))]
lta_mapping_lta2019_onetomany_ids <- lta_mapping %>% filter(lta2019_id %in% lta2019_onetomany_ids) %>% filter(!duplicated(lta2019_id))
lta_mapping_lta2019_onetoone_ids <- lta_mapping %>% filter(!lta2019_id %in% lta2019_onetomany_ids)
lta_mapping <- rbind(lta_mapping_lta2019_onetoone_ids, lta_mapping_lta2019_onetomany_ids)

# Function to extract all clades from a tree
get_clades <- function(tree) {
  bipart <- prop.part(tree)
  clades <- lapply(bipart, function(part) {
    return(sort(tree$tip.label[part]))
  })
  return(clades)
}

# Function to find parent proteins for a given set of child proteins (by overlap)
find_parent <- function(query_protein_ids, target_protein_ids) {
  query_protein_ids <- unique(unlist(strsplit(query_protein_ids, split=",")))
  shared_protein_ids <- intersect(target_protein_ids, query_protein_ids)
  return(length(shared_protein_ids) > 0)
}


for (k in 1:length(OG_list)) {
  selected_OG_id <- OG_list[k]
  if (bool_verbose) {
    print(paste0(selected_OG_id, " (", k, "/", length(OG_list), ")"))
  }
  
  outtablefile <- paste0(out_dir, "/", selected_tax_level, "/tsv/", selected_OG_id, "_consensus50_blopt_normalized_branch_length.tsv")
  if (bool_skip_completed & file.exists(outtablefile)) {
    print("Already completed.")
    next
  }
  
  # Read in unrooted branch-length optimized consensus tree
  consensus_blopt_tree_filename <- paste0(base_dir, "/", selected_OG_id, "_consensus50_blopt.treefile")
  consensus_blopt_iqtree_filename <- paste0(base_dir, "/", selected_OG_id, "_consensus50_blopt.iqtree")
  if (!file.exists(consensus_blopt_tree_filename) | !file.exists(consensus_blopt_iqtree_filename)) {
    print(paste0("Tree not completed: ", selected_OG_id, "_consensus50_blopt.treefile"))
    next
  }
  consensus_blopt_tree <- read.tree(consensus_blopt_tree_filename)
  
  # Map new Lta accessions and keep only the identical proteins in ltaref (5689)
  all_protein_ids <- consensus_blopt_tree$tip.label
  lta2019_protein_ids <- all_protein_ids[grepl("^5689_GET", all_protein_ids)]
  if (length(lta2019_protein_ids) > 0) {
    lta2019_protein_ids_remove <- lta2019_protein_ids[!lta2019_protein_ids %in% lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)]]
    lta2019_protein_ids_retain <- lta2019_protein_ids[lta2019_protein_ids %in% lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)]]
    
    # Remove lta2019 proteins that are absent in ltaref
    if (length(lta2019_protein_ids_remove) > 0) {
      consensus_blopt_tree <- drop.tip(consensus_blopt_tree, lta2019_protein_ids_remove)
    }
    
    # Replace lta2019 protein ids with ltaref protein ids
    if (length(lta2019_protein_ids_retain) > 0) {
      consensus_blopt_tree$tip.label[which(consensus_blopt_tree$tip.label %in% lta2019_protein_ids_retain)] <- lta_mapping$ltaref_id[match(consensus_blopt_tree$tip.label[which(consensus_blopt_tree$tip.label %in% lta2019_protein_ids_retain)], lta_mapping$lta2019_id)]
    }
  }
  
  # Get parent and child PhROGs@LECA 
  phrogs_long_eukaryota_parent_curr <- phrogs_long_eukaryota_parent %>% filter(OG_id == selected_OG_id)
  phrogs_long_eukaryota_curr <- phrogs_long_eukaryota %>% filter(OG_id == selected_OG_id)
  # Remove any proteins not found in tree (dropped Lta proteins)
  phrogs_long_eukaryota_parent_curr <- phrogs_long_eukaryota_parent_curr %>% filter(protein_id %in% consensus_blopt_tree$tip.label)
  phrogs_long_eukaryota_curr <- phrogs_long_eukaryota_curr %>% filter(protein_id %in% consensus_blopt_tree$tip.label)
  
  # Exclude nonvertical proteins
  if (bool_exclude_nonvertical) {
    phrogs_long_eukaryota_parent_curr_nonvertical <- phrogs_long_eukaryota_parent_curr %>% filter(BOOL_NONVERTICAL)
    phrogs_long_eukaryota_curr_nonvertical <- phrogs_long_eukaryota_curr %>% filter(BOOL_NONVERTICAL)
    nonvertical_protein_ids <- c(phrogs_long_eukaryota_parent_curr_nonvertical$protein_id, phrogs_long_eukaryota_curr_nonvertical$protein_id)
    
    phrogs_long_eukaryota_parent_curr <- phrogs_long_eukaryota_parent_curr %>% filter(!protein_id %in% nonvertical_protein_ids)
    phrogs_long_eukaryota_curr <- phrogs_long_eukaryota_curr %>% filter(!protein_id %in% nonvertical_protein_ids)
    if (length(nonvertical_protein_ids) > 0) {
      consensus_blopt_tree <- drop.tip(consensus_blopt_tree, nonvertical_protein_ids)
    }
  }
  
  # Prune small eukaryotic clades
  if (bool_exclude_small_clades) {
    phrogs_long_eukaryota_parent_curr_small_clades <- phrogs_long_eukaryota_parent_curr %>% group_by(PROG_id) %>% filter(length(unique(taxid)) < n_minimum_species_leca)
    phrogs_long_eukaryota_curr_small_clades <- phrogs_long_eukaryota_curr %>% group_by(PROG_id) %>% filter(length(unique(taxid)) < n_minimum_species_leca)
    small_clade_protein_ids <- unique(c(phrogs_long_eukaryota_parent_curr_small_clades$protein_id, phrogs_long_eukaryota_curr_small_clades$protein_id))
    
    phrogs_long_eukaryota_parent_curr <- phrogs_long_eukaryota_parent_curr %>% filter(!protein_id %in% small_clade_protein_ids)
    phrogs_long_eukaryota_curr <- phrogs_long_eukaryota_curr %>% filter(!protein_id %in% small_clade_protein_ids)
    if (length(small_clade_protein_ids) > 0) {
      if (length(small_clade_protein_ids) == Ntip(consensus_blopt_tree)) {
        print(paste0("Quitting since no clades pass minimum species threshold for ", selected_OG_id))
        quit(save = "no")
      } else {
        consensus_blopt_tree <- drop.tip(consensus_blopt_tree, small_clade_protein_ids)
      }
    }
  }
  
  phrogs_long_eukaryota_parent_leca <- phrogs_long_eukaryota_parent_curr %>% filter(label %in% ancestral_tax_labels)
  
  phrogs_long_eukaryota_curr_leca_wide <- phrogs_long_eukaryota_curr %>% filter(label == selected_tax_level) %>% group_by(PROG_id, OG_id) %>% summarize(protein_ids = paste0(protein_id, collapse=","), .groups = "drop")
  phrogs_long_eukaryota_curr_all_wide <- phrogs_long_eukaryota_curr %>% group_by(PROG_id, OG_id, label) %>% summarize(protein_ids = paste0(protein_id, collapse=","), .groups = "drop")
  
  if (nrow(phrogs_long_eukaryota_curr_leca_wide) == 0) {
    print(paste0("Warning: No LECA PhROGs found for ", selected_OG_id))
    next
  }
  
  ## Get stem length and median euk length for all PhROGs. Iterate through each parent PhROG and all their children
  PROG_list <- unique(phrogs_long_eukaryota_parent_leca$PROG_id)
  bl_df <- c()
  for (i in 1:length(PROG_list)) {
    parent_PROG_id <- PROG_list[i]
    if (bool_verbose) {
      print(paste0(parent_PROG_id, " (", i, "/", length(PROG_list), ")"))
    }
    
    # Find LECA child clades
    phrogs_long_eukaryota_parent_leca_curr <- phrogs_long_eukaryota_parent_leca %>% filter(PROG_id == parent_PROG_id)
    curr_parent_clade_protein_ids <- phrogs_long_eukaryota_parent_leca_curr$protein_id
    leca_child_indexes <- which(sapply(phrogs_long_eukaryota_curr_leca_wide$protein_ids, find_parent, curr_parent_clade_protein_ids))
    if (length(leca_child_indexes) == 0) {
      print(paste0("Warning: No child clade found for ", parent_PROG_id))
      next
    }
    # Find all child clades
    all_child_indexes <- which(sapply(phrogs_long_eukaryota_curr_all_wide$protein_ids, find_parent, curr_parent_clade_protein_ids))
    
    # Get full set of proteins (including non-vertical) for rooting
    phrogs_long_eukaryota_parent_curr <- phrogs_long_eukaryota_parent %>% filter(label %in% ancestral_tax_labels) %>% filter(OG_id == selected_OG_id) %>% filter(PROG_id == parent_PROG_id)
    curr_parent_clade_all_protein_ids <- phrogs_long_eukaryota_parent_curr$protein_id
    
    ## Root consensus tree 
    mrca_index <- get_mrca_of_set(consensus_blopt_tree, curr_parent_clade_all_protein_ids)
    outgroup_protein_ids <- consensus_blopt_tree$tip.label[!consensus_blopt_tree$tip.label %in% curr_parent_clade_all_protein_ids]
    if (length(outgroup_protein_ids) > 0) {
      # Root tree at deepest outgroup protein
      outgroup_protein_edge_distances <- get_pairwise_distances(consensus_blopt_tree, outgroup_protein_ids, rep(mrca_index, length(outgroup_protein_ids)), as_edge_counts = TRUE)
      outgroup_max_distance_tiplabels <- outgroup_protein_ids[outgroup_protein_edge_distances == max(outgroup_protein_edge_distances)]
      
      if (bool_verbose) {
        print(paste0("Rooting at deepest protein in outgroup: ", outgroup_max_distance_tiplabels[1], sep=","))
      }
      consensus_blopt_tree_claderoot <- root(consensus_blopt_tree, outgroup = outgroup_max_distance_tiplabels[1], resolve.root = TRUE)
      
    } else {
      # If no outgroup proteins available, use root from raw consensus tree
      
      # Read in raw alerax consensus tree (rooted)
      raw_consensus_tree_filename <- paste0(alerax_dir, "/", selected_OG_id, "/reconciliations/summaries/family_1_consensus_50.newick")
      raw_consensus_tree <- read.tree(raw_consensus_tree_filename)
      
      # Get the monophyletic clades descending from the root of raw consensus tree
      root_daughter_nodes <- raw_consensus_tree$edge[which(raw_consensus_tree$edge[,1] == find_root(raw_consensus_tree)), 2]
      if (length(root_daughter_nodes) == 2) {
        # Raw consensus tree has binary root
        if (root_daughter_nodes[1] > Ntip(raw_consensus_tree)) {
          raw_consensus_tree_rooted_clade_subtree <- get_subtree_at_node(raw_consensus_tree, root_daughter_nodes[1] - Ntip(raw_consensus_tree))$subtree
          root_monophyletic_clade_tiplabels_daughter1 <- raw_consensus_tree_rooted_clade_subtree$tip.label
        } else {
          root_monophyletic_clade_tiplabels_daughter1 <- raw_consensus_tree$tip.label[root_daughter_nodes[1]]
        }
        if (root_daughter_nodes[2] > Ntip(raw_consensus_tree)) {
          raw_consensus_tree_rooted_clade_subtree <- get_subtree_at_node(raw_consensus_tree, root_daughter_nodes[2] - Ntip(raw_consensus_tree))$subtree
          root_monophyletic_clade_tiplabels_daughter2 <- raw_consensus_tree_rooted_clade_subtree$tip.label
        } else {
          root_monophyletic_clade_tiplabels_daughter2 <- raw_consensus_tree$tip.label[root_daughter_nodes[2]]
        }
        # Take the smaller root clade as outgroup
        if (length(root_monophyletic_clade_tiplabels_daughter1) <= length(root_monophyletic_clade_tiplabels_daughter2)) {
          root_monophyletic_clade_tiplabels <- root_monophyletic_clade_tiplabels_daughter1
          nonroot_monophyletic_clade_tiplabels <- root_monophyletic_clade_tiplabels_daughter2
        } else {
          root_monophyletic_clade_tiplabels <- root_monophyletic_clade_tiplabels_daughter2
          nonroot_monophyletic_clade_tiplabels <- root_monophyletic_clade_tiplabels_daughter1
        }
      } else {
        # Raw consensus tree has multifurcating root. Determine whether LECA clades are within a resolved bipartition of the raw consensus tree
        leca_clades_tiplabels <- unlist(strsplit(phrogs_long_eukaryota_curr_leca_wide$protein_ids[leca_child_indexes], split=","))
        
        mrca_index <- get_mrca_of_set(raw_consensus_tree, leca_clades_tiplabels)
        if (mrca_index == find_root(raw_consensus_tree)) {
          print(paste0("Warning: Uncertain root position in raw consensus tree with respect to LECA clades. ", parent_PROG_id))
          # Root on largest clade within raw consensus tree
          raw_consensus_tree_clades <- get_clades(raw_consensus_tree)
          raw_consensus_tree_clades_count <- unlist(lapply(raw_consensus_tree_clades, function(lst) length(unique(lst))))
          largest_subclade_index <- order(raw_consensus_tree_clades_count, decreasing=TRUE)[2]
          root_monophyletic_clade_tiplabels <- unlist(raw_consensus_tree_clades[largest_subclade_index])
          nonroot_monophyletic_clade_tiplabels <- raw_consensus_tree$tip.label[!raw_consensus_tree$tip.label %in% root_monophyletic_clade_tiplabels]
          # next
        } else {
          print(paste0("Warning: Raw consensus tree has multifurcating root. ", parent_PROG_id))
          root_monophyletic_clade_tiplabels <- get_subtree_at_node(raw_consensus_tree, mrca_index - Ntip(raw_consensus_tree))$subtree$tip.label
          nonroot_monophyletic_clade_tiplabels <- raw_consensus_tree$tip.label[!raw_consensus_tree$tip.label %in% root_monophyletic_clade_tiplabels]
        }
      }
      
      # Map outgroup clade tips to the new tree for rooting
      consensus_blopt_tree_claderoot <- consensus_blopt_tree
      new_root_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, root_monophyletic_clade_tiplabels)
      # If root is at the original default root, cannot directly root due to ambiguity when resolving a node in an unrooted tree with no explicit outgroup
      if (new_root_index == find_root(consensus_blopt_tree_claderoot)) {
        if (bool_verbose) {
          print("New root same as default root. Rerooting...")
        }
        
        # First reroot at a nonroot monophyletic clade, then reroot
        temp_root_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, nonroot_monophyletic_clade_tiplabels)
        if (temp_root_index > Ntip(consensus_blopt_tree_claderoot)) {
          # Root is at a node
          consensus_blopt_tree_claderoot <- root(consensus_blopt_tree_claderoot, node = temp_root_index, resolve.root = TRUE)
        } else {
          # Root is at a tip
          consensus_blopt_tree_claderoot <- root(consensus_blopt_tree_claderoot, outgroup = consensus_blopt_tree_claderoot$tip.label[temp_root_index], resolve.root = TRUE)
        }
        
        new_root_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, root_monophyletic_clade_tiplabels)
      }
      
      # Root the tree
      if (new_root_index > Ntip(consensus_blopt_tree_claderoot)) {
        # Root is at a node
        if (bool_verbose) {
          print("Rooting on branch...")
        }
        
        # Root using phytools:reroot to root at the midpoint of the branch
        root_edge_length <- consensus_blopt_tree_claderoot$edge.length[consensus_blopt_tree_claderoot$edge[,2] == new_root_index]
        consensus_blopt_tree_claderoot <- reroot(consensus_blopt_tree_claderoot, new_root_index, position = root_edge_length/2)
        
      } else {
        # Root is at a tip
        if (bool_verbose) {
          print(paste0("Rooting on tip: ", consensus_blopt_tree_claderoot$tip.label[new_root_index]))
        }
        # consensus_blopt_tree_claderoot <- root(consensus_blopt_tree_claderoot, outgroup = new_root_index, resolve.root = TRUE)
        root_edge_length <- consensus_blopt_tree_claderoot$edge.length[consensus_blopt_tree_claderoot$edge[,2] == new_root_index]
        consensus_blopt_tree_claderoot <- reroot(consensus_blopt_tree_claderoot, new_root_index, position = root_edge_length/2)
      }
      
    }
    
    # Remove node support labels
    consensus_blopt_tree_claderoot$node.label <- rep("", length(consensus_blopt_tree_claderoot$node.label))
    
    # Get the acquisition node for the current clade
    parent_mrca_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, curr_parent_clade_protein_ids)
    acquisition_node_index <- consensus_blopt_tree_claderoot$edge[which(consensus_blopt_tree_claderoot$edge[,2] == parent_mrca_index), 1]
    
    # Get all clade node MRCAs
    clade_mrca_indexes <- c()
    for (j in 1:length(all_child_indexes)) {
      child_index <- all_child_indexes[j]
      phrogs_long_eukaryota_leca_curr_wide_child <- phrogs_long_eukaryota_curr_all_wide[child_index,]
      
      curr_clade_protein_ids <- unlist(strsplit(phrogs_long_eukaryota_leca_curr_wide_child$protein_ids, split=","))
      clade_node_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, curr_clade_protein_ids)
      clade_mrca_indexes <- c(clade_mrca_indexes, clade_node_index)
      
      if (clade_node_index - Ntip(consensus_blopt_tree_claderoot) > 0) {
        consensus_blopt_tree_claderoot$node.label[clade_node_index - Ntip(consensus_blopt_tree_claderoot)] <- paste0(phrogs_long_eukaryota_leca_curr_wide_child$label, "_", gsub(".*_", "", phrogs_long_eukaryota_leca_curr_wide_child$PROG_id))
      }
    }
    
    for (j in 1:length(leca_child_indexes)) {
      child_index <- leca_child_indexes[j]
      phrogs_long_eukaryota_leca_curr_wide_child <- phrogs_long_eukaryota_curr_leca_wide[child_index,]
      child_PROG_id <- phrogs_long_eukaryota_leca_curr_wide_child$PROG_id
      
      if (bool_verbose) {
        print(paste0("Child ", child_PROG_id, " (", j, "/", length(leca_child_indexes), ")"))
      }
      
      curr_clade_protein_ids <- unlist(strsplit(phrogs_long_eukaryota_leca_curr_wide_child$protein_ids, split=","))
      leca_node_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, curr_clade_protein_ids)
      
      # Get median eukaryote branch length
      if (bool_verbose) {
        print("Getting median eukaryote branch length...")
      }
      euk_clade_bl <- get_pairwise_distances(consensus_blopt_tree_claderoot, curr_clade_protein_ids, rep(leca_node_index, length(curr_clade_protein_ids)))
      if (bool_verbose) {
        print(paste0("Median euk clade bl: ", median(euk_clade_bl)))
        print(paste0("Mean euk clade bl: ", mean(euk_clade_bl)))
        print(paste0("SD euk clade bl: ", sd(euk_clade_bl)))
      }
      
      # Get stem length
      stem_bl <- NA
      # Check that acquisition node exists. May not exist if all tree proteins are in parent PhROG
      if (length(acquisition_node_index) == 1) {
        if (bool_verbose) {
          print("Getting stem length...")
        }
        stem_bl <- get_pairwise_distances(consensus_blopt_tree_claderoot, acquisition_node_index, leca_node_index)
        if (bool_verbose) {
          print(paste0("Stem length:", stem_bl))
          print(paste0("Normalized stem length:", stem_bl / median(euk_clade_bl)))
        }
        consensus_blopt_tree_claderoot$node.label[acquisition_node_index - Ntip(consensus_blopt_tree_claderoot)] <- parent_PROG_id
      }
      
      consensus_blopt_tree_claderoot$node.label[leca_node_index - Ntip(consensus_blopt_tree_claderoot)] <- paste0(selected_tax_level, "_", gsub(".*_", "", child_PROG_id))
      
      # Get duplication length
      duplication_bl <- NA
      oldest_duplication_bl <- NA
      if (length(all_child_indexes) >=2) {
        # Get duplication node: MRCA between current clade and any other LECA clade
        other_leca_node_indexes <- clade_mrca_indexes[clade_mrca_indexes != leca_node_index]
        candidate_duplication_mrcas <- get_pairwise_mrcas(consensus_blopt_tree_claderoot, rep(leca_node_index, length(other_leca_node_indexes)), other_leca_node_indexes)
        candidate_duplication_mrcas <- unique(candidate_duplication_mrcas)
        candidate_duplication_mrcas_distances <- get_pairwise_distances(consensus_blopt_tree_claderoot, rep(leca_node_index, length(candidate_duplication_mrcas)), candidate_duplication_mrcas, as_edge_counts = TRUE)
        
        # Assign the closest MRCA node as the duplication node 
        duplication_node_index <- candidate_duplication_mrcas[candidate_duplication_mrcas_distances == min(candidate_duplication_mrcas_distances)]
        
        if (length(duplication_node_index) == 1) {
          duplication_bl <- get_pairwise_distances(consensus_blopt_tree_claderoot, duplication_node_index, leca_node_index)
          if (bool_verbose) {
            print(paste0("Duplication length:", duplication_bl))
            print(paste0("Normalized duplication length:", duplication_bl / median(euk_clade_bl)))
          }
          
          consensus_blopt_tree_claderoot$node.label[duplication_node_index - Ntip(consensus_blopt_tree_claderoot)] <- "Duplication"
        }
        
        # Get the oldest duplication node 
        oldest_duplication_node_index <- candidate_duplication_mrcas[candidate_duplication_mrcas_distances == max(candidate_duplication_mrcas_distances)]
        if (length(oldest_duplication_node_index) == 1) {
          oldest_duplication_bl <- get_pairwise_distances(consensus_blopt_tree_claderoot, oldest_duplication_node_index, leca_node_index)
        }
      }
      
      curr_bl_df <- data.frame(OG_id = selected_OG_id, parent_PhROG_id = parent_PROG_id, child_PhROG_id = child_PROG_id, median_euk_clade_bl = median(euk_clade_bl), mean_euk_clade_bl = mean(euk_clade_bl), sd_euk_clade_bl = sd(euk_clade_bl), stem_bl, normalized_stem_bl = stem_bl / median(euk_clade_bl), duplication_bl, normalized_duplication_bl = duplication_bl / median(euk_clade_bl), oldest_duplication_bl, normalized_oldest_duplication_bl = oldest_duplication_bl / median(euk_clade_bl))
      bl_df <- rbind(bl_df, curr_bl_df)
    }
    
    # # Write out tree with node labels
    # write.tree(consensus_blopt_tree_claderoot, here("reconciled_consensus_trees_for_LECA_timing", parent_PROG_id, "_consensus50_blopt.treefile.rooted"))
  }
  
  # # Write out branch lengths table
  # write.table(bl_df, here("data/branch_length_timing", paste0("branch_lengths_", selected_tax_level, "_combined.tsv")), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
  
}

