# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(castor))
suppressMessages(library(phytools))

# Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
phrog_base_dir <- args[1]
consensus50_blopt_tree_dir <- args[2]
alerax_dir <- args[3]
out_dir <- args[4]
OG_id <- args[5]
selected_tax_levels_file <- args[6]
bool_verbose <- FALSE


### Initialize parameters and data
# Set random seed
set.seed(42)

# Focal taxonomic levels for which to generate PhROGs
selected_tax_levels <- read.table(here("data/phylogenetically_resolved_orthogroups", "selected_tax_levels_for_phrogs_wholeproteome.txt"), header=FALSE)$V1

# Core species taxids
completed_mitoproteomes_species_list <- c("9606", "559292", "3702", "1257118", "5689", "185431", "5741", "32595")

consensus50_bl_opt_tree_filename <- here("reconciled_consensus_trees_with_supports", paste0(OG_id, "_consensus50_blopt.treefile"))
if (!file.exists(consensus50_bl_opt_tree_filename)) {
  print(paste0("Tree file not found: ", consensus50_bl_opt_tree_filename))
  quit(save="no")
} else {
  consensus_blopt_tree <- read.tree(consensus50_bl_opt_tree_filename)
}

# Read in Eukaryota parent PhROGs for the current OG
parent_progs_long_raw <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(parent_progs_long_raw) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_parent <- parent_progs_long_raw %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% filter(OG_id == OG_id) %>% filter(!duplicated(protein_id))

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


### Determine ancestral localizations
# Read in experimental and mtDNA mito proteins
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_2025.09.30.tsv"), sep="\t", header=TRUE)

# Read in deeploc predictions
deeploc_results <- read.table(file.path(supplemental_data_directory, "TableS7_retrained_DeepLoc_predictions.tsv"), header=TRUE)
colnames(deeploc_results) <- c("Protein_ID", "Mitochondrion")

# Read in organelle-encoded proteins
mtdna_proteins <- read.table(here("data/deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1
nonmito_organelle_proteins <- read.table(here("data/deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
missing_mtdna_taxids <- read.table(here("data/deeploc", "missing_mtdna_taxids.txt"))$V1
missing_nonmito_organelle_taxids <- read.table(here("data/deeploc", "missing_nonmito_organelle_taxids.txt"))$V1

# Function to extract all clades from a tree
get_clades <- function(tree) {
  bipart <- prop.part(tree)
  clades <- lapply(bipart, function(part) {
    return(sort(tree$tip.label[part]))
  })
  return(clades)
}

## Iterate through parent PhROGs
mito_likelihood_tips_agg <- c()
parent_phrogs <- unique(phrogs_long_eukaryota_parent$PROG_id)
for (i in 1:length(parent_phrogs)) {
  curr_parent_phrog <- parent_phrogs[i]
  
  if (bool_verbose) {
    print(paste0(curr_parent_phrog, ": ", i, " of ", length(parent_phrogs)))
  }
  
  phrogs_long_eukaryota_parent_curr <- phrogs_long_eukaryota_parent %>% filter(PROG_id == curr_parent_phrog)
  
  # Get full set of proteins (including non-vertical) for rooting
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
    raw_consensus_tree_filename <- paste0(alerax_dir, "/", OG_id, "/reconciliations/summaries/family_1_consensus_50.newick")
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
      # Raw consensus tree has multifurcating root. Determine whether parent Eukaryota clades are within a resolved bipartition of the raw consensus tree
      mrca_index <- get_mrca_of_set(raw_consensus_tree, curr_parent_clade_all_protein_ids)
      if (mrca_index == find_root(raw_consensus_tree)) {
        print(paste0("Warning: Uncertain root position in raw consensus tree with respect to parent Eukaryota clades. ", curr_parent_phrog))
        # Root on largest clade within raw consensus tree
        raw_consensus_tree_clades <- get_clades(raw_consensus_tree)
        raw_consensus_tree_clades_count <- unlist(lapply(raw_consensus_tree_clades, function(lst) length(unique(lst))))
        largest_subclade_index <- order(raw_consensus_tree_clades_count, decreasing=TRUE)[2]
        root_monophyletic_clade_tiplabels <- unlist(raw_consensus_tree_clades[largest_subclade_index])
        nonroot_monophyletic_clade_tiplabels <- raw_consensus_tree$tip.label[!raw_consensus_tree$tip.label %in% root_monophyletic_clade_tiplabels]
        # next
      } else {
        print(paste0("Warning: Raw consensus tree has multifurcating root. ", curr_parent_phrog))
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
      
      # First reroot at a random tip in clade, then reroot
      temp_root_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, nonroot_monophyletic_clade_tiplabels)
      consensus_blopt_tree_claderoot <- root(consensus_blopt_tree_claderoot, node = temp_root_index, resolve.root = TRUE)
      
      new_root_index <- get_mrca_of_set(consensus_blopt_tree_claderoot, root_monophyletic_clade_tiplabels)
    }
    
    # Root the tree
    if (new_root_index > Ntip(consensus_blopt_tree_claderoot)) {
      # Root is at a node
      if (bool_verbose) {
        print("Rooting on branch...")
      }
      # consensus_blopt_tree_claderoot <- root(consensus_blopt_tree_claderoot, node = new_root_index, resolve.root = TRUE)
      
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
  
  curr_tree <- keep.tip(consensus_blopt_tree_claderoot, curr_parent_clade_all_protein_ids[curr_parent_clade_all_protein_ids %in% consensus_blopt_tree_claderoot$tip.label])
  
  # Initialize
  protein_accessions <- curr_tree$tip.label
  curr_tree_mito_localization_states <- rep(NA, Ntip(curr_tree))
  
  # Mark states for reference species with complete localization data: binary mito or not mito
  curr_tree_mito_localization_states[which(gsub("_.*", "", curr_tree$tip.label) %in% completed_mitoproteomes_species_list)] <- 0 # nonmito
  curr_tree_mito_localization_states[which(protein_accessions %in% gold_gene_accession_OG_id_df$gene_accession[gold_gene_accession_OG_id_df$taxid %in% completed_mitoproteomes_species_list])] <- 1 # mito
  
  ## Mk model on binary label.
  # Initialize priors matrix: column 1 is nonmito probability, column 2 is mito probability. Each row sums to 1
  default_prior <- NA # no prior information for unannotated species
  tip_priors_mat <- matrix(data = default_prior, nrow = length(curr_tree_mito_localization_states), ncol = 2)
  
  # Use deeploc predictions as priors
  tip_priors_mat[which(protein_accessions %in% deeploc_results$Protein_ID),2] <- deeploc_results$Mitochondrion[match(protein_accessions[which(protein_accessions %in% deeploc_results$Protein_ID)], deeploc_results$Protein_ID)]
  
  # Assign binary labels to organelle-encoded proteins.
  tip_priors_mat[which(protein_accessions %in% nonmito_organelle_proteins),2] <- 0
  tip_priors_mat[which(protein_accessions %in% mtdna_proteins),2] <- 1
  curr_tree_mito_localization_states[which(protein_accessions %in% nonmito_organelle_proteins)] <- 0
  curr_tree_mito_localization_states[which(protein_accessions %in% mtdna_proteins)] <- 1
  
  # If there are any organelle-encoded proteins present, omit deeploc predictions for species without mtDNA annotation since we don't know which proteins are nuclear-encoded versus mtDNA-encoded
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
  
  # Relabel priors for the MitoCarta species proteins with binary labels.
  tip_priors_mat[which(curr_tree_mito_localization_states == 0),2] <- 0
  tip_priors_mat[which(curr_tree_mito_localization_states == 1),2] <- 1
  
  # Ignore unedited mtDNA proteins and likely NUMTs in Arabidopsis
  ath_proteins_curr <- grep("3702_", curr_tree$tip.label, fixed=TRUE, value=TRUE)
  ath_uneditedmtdna_or_numt_proteins <- read.table(here("data/mito_orthogroups", "ath_unedited_numt_mtdna_proteins.txt"))$V1
  ath_uneditedmtdna_or_numt_proteins <- ath_uneditedmtdna_or_numt_proteins[ath_uneditedmtdna_or_numt_proteins %in% ath_proteins_curr]
  if (length(ath_uneditedmtdna_or_numt_proteins) > 0) {
    if (bool_verbose) {
      print(paste0("Ignoring Arabidopsis unedited mtDNA proteins: ", paste0(ath_uneditedmtdna_or_numt_proteins, collapse=",")))
    }
    
    curr_tree_mito_localization_states[which(protein_accessions %in% ath_uneditedmtdna_or_numt_proteins)] <- NA
    tip_priors_mat[which(protein_accessions %in% ath_uneditedmtdna_or_numt_proteins),2] <- NA
  }
  
  # Mark nonmito state as the complement of mito state
  tip_priors_mat[,1] <- 1 - tip_priors_mat[,2]
  
  if (sum(!is.na(tip_priors_mat[,2])) > 1) {
    # If there are multiple labels available
    mito_likelihood_mk <- rep(NA, Ntip(curr_tree) + Nnode(curr_tree))
    
    # Mk model with equal rates. Use trycatch to handle model fit failure
    tryCatch(
      {
        result <- hsp_mk_model(tree=curr_tree, tip_states=NULL, tip_priors=tip_priors_mat, Nstates=2, rate_model = "ER", Ntrials=5)
        mito_likelihood_mk <- result$likelihoods[,2] # take the mito state likelihood
      },
      error = function(cond) {
        # Code to execute if an error occurs
        print(paste0("Error in hsp_mk_model: ", curr_parent_phrog))
        mito_likelihood_mk <- rep(NA, Ntip(curr_tree) + Nnode(curr_tree))
      },
      warning = function(cond) {
        # Code to execute if a warning occurs
        print(paste0("Warning in hsp_mk_model: ", curr_parent_phrog))
        mito_likelihood_mk <- rep(NA, Ntip(curr_tree) + Nnode(curr_tree))
      }
    )
    
  } else if (sum(!is.na(tip_priors_mat[,2])) == 1) {
    # If there is only 1 label available, just use its score to inpute the rest
    mito_likelihood_mk <- rep(tip_priors_mat[!is.na(tip_priors_mat[,2]),2], Ntip(curr_tree) + Nnode(curr_tree))
  } else {
    # If there are no scores available, assume nonmito and set to 0
    print(paste0("No annotated species proteins present in ", OG_id))
    mito_likelihood_mk <- rep(0, Ntip(curr_tree) + Nnode(curr_tree))
  }
  
  ## Parsimony with equal transition cost only on certain labels (MitoCarta proteins + organelle-encoded)
  cost_matrix <- matrix(data=c(0,1,1,0), nrow=2, ncol=2)
  if (sum(!is.na(curr_tree_mito_localization_states)) > 1) {
    # If there are multiple labels available
    result <- hsp_max_parsimony(curr_tree, curr_tree_mito_localization_states+1, Nstates=2, transition_costs = cost_matrix, edge_exponent=0)
    mito_likelihood_parsimony <- result$likelihoods[,2] # take the mito likelihood
  } else if (sum(!is.na(curr_tree_mito_localization_states)) == 1) {
    # If there is only one label available, use it to impute the rest
    mito_likelihood_parsimony <- rep(curr_tree_mito_localization_states[!is.na(curr_tree_mito_localization_states)], Ntip(curr_tree) + Nnode(curr_tree))
  } else {
    # If there are no scores available, assume nonmito and set to 0
    mito_likelihood_parsimony <- rep(0, Ntip(curr_tree) + Nnode(curr_tree))
  }
  
  # Get mito localization predictions at tips
  mito_likelihood_tips <- data.frame(OG_id = OG_id, protein_id = curr_tree$tip.label, mito_localization_prob = tip_priors_mat[,2], mito_localization_prob_mk = mito_likelihood_mk[1:Ntip(curr_tree)], mito_localization_prob_parsimony = mito_likelihood_parsimony[1:Ntip(curr_tree)])
  mito_likelihood_tips_agg <- rbind(mito_likelihood_tips_agg, mito_likelihood_tips)
  
  # Add ancestral localizations
  curr_tree_mk <- curr_tree
  curr_tree_mk$node.label <- mito_likelihood_mk[(Ntip(curr_tree)+1):(Ntip(curr_tree)+Nnode(curr_tree))]
  curr_tree_parsimony <- curr_tree
  curr_tree_parsimony$node.label <- mito_likelihood_parsimony[(Ntip(curr_tree)+1):(Ntip(curr_tree)+Nnode(curr_tree))]
  
  if (i == 1) {
    trees_parent_phrogs_mk <- list(curr_tree_mk)
    trees_parent_phrogs_parsimony <- list(curr_tree_parsimony)
  } else {
    trees_parent_phrogs_mk <- c(trees_parent_phrogs_mk, list(curr_tree_mk))
    trees_parent_phrogs_parsimony <- c(trees_parent_phrogs_parsimony, list(curr_tree_parsimony))
  }
}

## Update PhROGs
for (i in 1:length(selected_tax_levels)) {
  selected_tax_level <- selected_tax_levels[i]
  if (bool_verbose) {
    print(paste0(selected_tax_level, ": ", i, " of ", length(selected_tax_levels)))
  }
  
  curr_out_dir <- paste0(out_dir, "/", selected_tax_level)
  dir.create(curr_out_dir, showWarnings = FALSE)
  outfilename <- paste0(curr_out_dir, "/", OG_id, "_", selected_tax_level, "_PhROGs_long.tsv")
  
  phrogs_long_filename <- paste0(phrog_base_dir, "/", selected_tax_level, "/", OG_id, "_", selected_tax_level, "_PhROGs_long.tsv")
  if (!file.exists(phrogs_long_filename)) {
    next
  } else if (file.size(phrogs_long_filename) == 0) {
    # No PhROGs at this taxonomic level. Create empty file
    file.create(outfilename)
    next
  } else {
    phrogs_long_original <- read.table(phrogs_long_filename, sep="\t", header=FALSE)
    colnames(phrogs_long_original) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
    phrogs_long_original$mito_localization_prob_mk <- NA
    phrogs_long_original$mito_localization_prob_parsimony <- NA
    
    phrogs_long <- phrogs_long_original %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) #%>% filter(!duplicated(protein_id))
  }
  
  PhROG_ids <- unique(phrogs_long$PROG_id)
  for (i in 1:length(PhROG_ids)) {
    PhROG_id_curr <- PhROG_ids[i]
    
    if (bool_verbose) {
      print(paste0(PhROG_id_curr, ": ", i, " of ", length(PhROG_ids)))
    }
    
    phrogs_long_curr <- phrogs_long %>% filter(PROG_id == PhROG_id_curr)
    
    matched_tree_index <- c()
    for (i in 1:length(trees_parent_phrogs_mk)) {
      trees_parent_phrogs_curr <- trees_parent_phrogs_mk[[i]]
      
      if (all(phrogs_long_curr$protein_id %in% trees_parent_phrogs_curr$tip.label)) {
        matched_tree_index <- c(matched_tree_index, i)
      }
    }
    
    if (length(matched_tree_index) == 1) {
      mrca_index <- get_mrca_of_set(trees_parent_phrogs_mk[[matched_tree_index]], phrogs_long_curr$protein_id)
      
      if (mrca_index <= Ntip(trees_parent_phrogs_mk[[matched_tree_index]])) {
        mito_likelihood_mk_curr <- mito_likelihood_tips_agg$mito_localization_prob_mk[match(phrogs_long_curr$protein_id, mito_likelihood_tips_agg$protein_id)]
        mito_likelihood_parsimony_curr <- mito_likelihood_tips_agg$mito_localization_prob_parsimony[match(phrogs_long_curr$protein_id, mito_likelihood_tips_agg$protein_id)]
      } else {
        mito_likelihood_mk_curr <- trees_parent_phrogs_mk[[matched_tree_index]]$node.label[mrca_index - Ntip(trees_parent_phrogs_mk[[matched_tree_index]])]
        mito_likelihood_parsimony_curr <- trees_parent_phrogs_parsimony[[matched_tree_index]]$node.label[mrca_index - Ntip(trees_parent_phrogs_parsimony[[matched_tree_index]])]
      }
      
      # Update localization probabilities
      phrogs_long_original$mito_localization_prob <- mito_likelihood_tips_agg$mito_localization_prob[match(phrogs_long_original$protein_id, mito_likelihood_tips_agg$protein_id)]
      phrogs_long_original$mito_localization_prob_mk[which(phrogs_long_original$PROG_id == PhROG_id_curr)] <- mito_likelihood_mk_curr
      phrogs_long_original$mito_localization_prob_parsimony[which(phrogs_long_original$PROG_id == PhROG_id_curr)] <- mito_likelihood_parsimony_curr
      
    } else if (length(matched_tree_index) == 0) {
      print(paste0("Warning: no matched tree index for ", PhROG_id_curr))
    } else if (length(matched_tree_index) > 1) {
      print(paste0("Warning: multiple matched tree indexes for ", PhROG_id_curr))
    }
  }
  
  # # Write out
  # write.table(phrogs_long_original, outfilename, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
