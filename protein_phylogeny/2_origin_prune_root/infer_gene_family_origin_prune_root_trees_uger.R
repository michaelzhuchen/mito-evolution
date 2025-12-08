### Infer domain of origin, prune prokaryotes, and root by MAD (in largest prokaryotic outgroup if prokaryotic origin)

# Load libraries. Also requires MAD-2.2
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(castor))
suppressMessages(library(stringr))

# Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
base_dir <- args[1]
out_dir <- args[2]
MSA_id <- args[3]
bool_verbose <- FALSE

# Set random seed
set.seed(42)

# Set thresholds
clade_purity_threshold <- 0.9
min_frac_prok_species_per_subgroup_for_prok_origin <- 0.25
min_frac_subgroups_per_prokgroup_for_prok_origin <- 0.5
min_frac_euk_species_for_euk_origin <- 0.9
min_prok_species_in_clade_to_not_prune <- 10

if (bool_verbose) {
  print("Reading in data")
}

# Identify ML tree: either treefile or contree
file_lines <- readLines(paste0(base_dir, "/", MSA_id, ".iqtree"))
bool_use_consensus_tree <- any(grepl("Consensus tree has higher likelihood than ML tree found", file_lines))
if (bool_use_consensus_tree) {
  ML_tree_suffix <- ".contree"
} else {
  ML_tree_suffix <- ".treefile"
}
ML_tree_filename <- paste0(base_dir, "/", MSA_id, ML_tree_suffix)

# Read in ufboot trees
ufboot_trees_filename <- paste0(base_dir, "/", MSA_id, ".ufboot")

if (!file.exists(ML_tree_filename) | !file.exists(ufboot_trees_filename)) {
  quit(save="no")
}

ML_tree <- read.tree(ML_tree_filename)
ufboot_trees <- read.tree(ufboot_trees_filename)
n_ufboot_trees <- length(ufboot_trees)

# Drop duplicated proteins. Muscle appends " dupelabel.X" to duplicated accessions
dup_proteins <- grep(" dupelabel", ML_tree$tip.label, fixed=TRUE, value=TRUE)
if (length(dup_proteins) > 0) {
  if (bool_verbose) {
    print(paste0("Dropping duplicated proteins: ", paste0(dup_proteins, collapse=",")))
  }
  ML_tree <- drop.tip(ML_tree, dup_proteins)
  for (i in 1:length(ufboot_trees)) {
    ufboot_trees[[i]] <- drop.tip(ufboot_trees[[i]], dup_proteins)
  }
}

# Drop mtDNA proteins from taxid 2086695 (incorrectly mapped sequences)
remove_mtDNA_proteins <- grep("2086695_", ML_tree$tip.label, fixed=TRUE, value=TRUE)
if (length(remove_mtDNA_proteins) > 0) {
  if (bool_verbose) {
    print(paste0("Dropping mtDNA proteins: ", paste0(remove_mtDNA_proteins, collapse=",")))
  }
  ML_tree <- drop.tip(ML_tree, remove_mtDNA_proteins)
  for (i in 1:length(ufboot_trees)) {
    ufboot_trees[[i]] <- drop.tip(ufboot_trees[[i]], remove_mtDNA_proteins)
  }
}

# Rename incorrect mtDNA taxids
rename_mtDNA_proteins <- grep("^192875_", ML_tree$tip.label, value=TRUE)
if (length(rename_mtDNA_proteins) > 0) {
  if (bool_verbose) {
    print(paste0("Rename mtDNA proteins: ", paste0(rename_mtDNA_proteins, collapse=",")))
  }
  ML_tree$tip.label[match(rename_mtDNA_proteins, ML_tree$tip.label)] <- gsub("^192875_", "595528_", rename_mtDNA_proteins)
  for (i in 1:length(ufboot_trees)) {
    ufboot_trees[[i]]$tip.label[match(rename_mtDNA_proteins, ufboot_trees[[i]]$tip.label)] <- gsub("^192875_", "595528_", rename_mtDNA_proteins)
  }
}


# Read in species tree
species_tree <- here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample.contree")
species_tree_labels <- c(species_tree$tip.label, species_tree$node.label)

# Read in uniprot tax data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
uniprot_proteomes_euk_taxids <- uniprot_proteomes_tax$TaxId[uniprot_proteomes_tax$domain == "Eukaryota"]

# Read in uniprot all species tax data
uniprot_proteomes_all_tax <- read.table(here("data/taxonomy", "uniprot_bacteria.downsample.level6_tax_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
bacteria_taxids <- uniprot_proteomes_all_tax$TaxId[which(uniprot_proteomes_all_tax$domain == "Bacteria")]
archaea_taxids <- uniprot_proteomes_all_tax$TaxId[which(uniprot_proteomes_all_tax$domain == "Archaea")]
n_bacteria_species <- length(bacteria_taxids)
n_archaea_species <- length(archaea_taxids)

# Get total species counts per prokaryotic subgroup
uniprot_proteomes_bacteria_subgroup_counts <- uniprot_proteomes_all_tax %>% filter(domain == "Bacteria") %>% group_by(tree_id, level_5) %>% summarize(n_species_per_subgroup = n(), .groups="keep") %>% mutate(n_species_per_subgroup_threshold = n_species_per_subgroup*min_frac_prok_species_per_subgroup_for_prok_origin)
uniprot_proteomes_archaea_subgroup_counts <- uniprot_proteomes_all_tax %>% filter(domain == "Archaea") %>% group_by(tree_id, kingdom) %>% summarize(n_species_per_subgroup = n(), .groups="keep") %>% mutate(n_species_per_subgroup_threshold = n_species_per_subgroup*min_frac_prok_species_per_subgroup_for_prok_origin)
colnames(uniprot_proteomes_bacteria_subgroup_counts)[2] <- "subgroup"
colnames(uniprot_proteomes_archaea_subgroup_counts)[2] <- "subgroup"
uniprot_proteomes_prokaryote_subgroup_counts <- rbind(uniprot_proteomes_bacteria_subgroup_counts, uniprot_proteomes_archaea_subgroup_counts)
uniprot_proteomes_prokaryote_subgroup_counts_summary <- uniprot_proteomes_prokaryote_subgroup_counts %>% group_by(tree_id) %>% summarize(n_subgroups_per_prokgroup = n(), .groups = "keep") %>% mutate(n_subgroups_per_prokgroup_threshold = n_subgroups_per_prokgroup*min_frac_subgroups_per_prokgroup_for_prok_origin)

# Get species counts per prokaryotic subgroup for largest prokaryotic clade
count_large_prokgroups_in_clade <- function(prok_species_in_clade) {
  prok_species_largest_prok_clade <- unique(unlist(strsplit(prok_species_in_clade, split=",")))
  prok_largest_clade_tax <- uniprot_proteomes_all_tax[which(uniprot_proteomes_all_tax$TaxId %in% prok_species_largest_prok_clade),]
  prok_largest_clade_bacteria_subgroup_counts <- prok_largest_clade_tax %>% filter(domain == "Bacteria") %>% group_by(tree_id, level_5) %>% summarize(n_species_per_subgroup = n(), .groups="keep")
  prok_largest_clade_archaea_subgroup_counts <- prok_largest_clade_tax %>% filter(domain == "Archaea") %>% group_by(tree_id, kingdom) %>% summarize(n_species_per_subgroup = n(), .groups="keep")
  colnames(prok_largest_clade_bacteria_subgroup_counts)[2] <- "subgroup"
  colnames(prok_largest_clade_archaea_subgroup_counts)[2] <- "subgroup"
  prok_largest_clade_prokaryote_subgroup_counts <- rbind(prok_largest_clade_bacteria_subgroup_counts, prok_largest_clade_archaea_subgroup_counts)
  
  # Identify prokaryotic subgroups passing the species thresholds
  prok_largest_clade_prokaryote_subgroup_counts$n_species_per_subgroup_threshold <- uniprot_proteomes_prokaryote_subgroup_counts$n_species_per_subgroup_threshold[match(prok_largest_clade_prokaryote_subgroup_counts$subgroup, uniprot_proteomes_prokaryote_subgroup_counts$subgroup)]
  prok_largest_clade_prokaryote_subgroup_counts_filter <- prok_largest_clade_prokaryote_subgroup_counts[which(prok_largest_clade_prokaryote_subgroup_counts$n_species_per_subgroup >= prok_largest_clade_prokaryote_subgroup_counts$n_species_per_subgroup_threshold),]
  prok_largest_clade_prokaryote_prokgroup_subgroups <- prok_largest_clade_prokaryote_subgroup_counts_filter %>% group_by(tree_id) %>% summarize(n_subgroups_per_prokgroup=n())
  
  # Identify prokaryotic groups passing the subgroup thresholds
  prok_largest_clade_prokaryote_prokgroup_subgroups$n_subgroups_per_prokgroup_threshold <- uniprot_proteomes_prokaryote_subgroup_counts_summary$n_subgroups_per_prokgroup_threshold[match(prok_largest_clade_prokaryote_prokgroup_subgroups$tree_id, uniprot_proteomes_prokaryote_subgroup_counts_summary$tree_id)]
  prok_largest_clade_prokaryote_prokgroup_filter <- prok_largest_clade_prokaryote_prokgroup_subgroups[which(prok_largest_clade_prokaryote_prokgroup_subgroups$n_subgroups_per_prokgroup >= prok_largest_clade_prokaryote_prokgroup_subgroups$n_subgroups_per_prokgroup_threshold),]
  return(nrow(prok_largest_clade_prokaryote_prokgroup_filter))
}

# Read in the prokaryote mmseqs2 cluster data
prokaryote_mmseqs2_clusters_taxids <- read.table(here("data/downsample_prokaryotes", "prokaryote_mmseqs2_clusters_taxids.tsv"), sep="\t", quote="", header=TRUE)
ML_proteins_remove_first_underscore <- sub("^[^_]*_", "", ML_tree$tip.label)
prokaryote_mmseqs2_clusters_taxids_in_tree <- prokaryote_mmseqs2_clusters_taxids[which(prokaryote_mmseqs2_clusters_taxids$rep_seq %in% ML_proteins_remove_first_underscore),]

if (bool_verbose) {
  print("Enumerating clades")
}

# Function to extract all clades from a tree
get_clades <- function(tree) {
  bipart <- prop.part(tree)
  clades <- lapply(bipart, function(part) {
    return(sort(tree$tip.label[part]))
  })
  return(clades)
}

# Iterate through all roots to find all possible clades
ML_tree_rooted_clades <- c()
for (node_num in 1:Nnode(ML_tree)) {
  # print(node_num)
  node_index <- node_num + Ntip(ML_tree)
  if (node_index != find_root(ML_tree)) {
    ML_tree_curr_rooted <- root(ML_tree, node = node_index, resolve.root = TRUE)
  } else {
    ML_tree_curr_rooted <- ML_tree
  }
  
  # ML_tree_rooted_subtrees <- subtrees(ML_tree_curr_rooted)
  ML_tree_rooted_clades_curr <- get_clades(ML_tree_curr_rooted)
  ML_tree_rooted_clades <- unique(c(ML_tree_rooted_clades, ML_tree_rooted_clades_curr))
}

# Add in single tips as possible clades
ML_tree_rooted_clades <- c(ML_tree_rooted_clades, ML_tree$tip.label)

## Find the largest monophyletic clades
if (bool_verbose) {
  print("Finding largest monophyletic clades")
}
monophyletic_clades_df <- c()

# Get eukaryote proteins and species counts
ML_tree_rooted_clades_species <- lapply(ML_tree_rooted_clades, function(lst) unique(gsub("_.*", "", lst)))
n_euk_species_in_clades <- unlist(lapply(ML_tree_rooted_clades_species, function(lst) sum(lst %in% uniprot_proteomes_euk_taxids)))
euk_species_in_clades <- unlist(lapply(ML_tree_rooted_clades_species, function(lst) paste0(lst[which(lst %in% uniprot_proteomes_euk_taxids)], collapse=",")))
euk_proteins_in_clades_raw <- lapply(ML_tree_rooted_clades, function(lst) lst[gsub("_.*", "", lst) %in% uniprot_proteomes_euk_taxids])
euk_proteins_in_clades_raw_collapse <- unlist(lapply(euk_proteins_in_clades_raw, function(lst) paste0(lst, collapse=",")))

# Get prokaryote proteins and species counts
prok_proteins_in_clades_raw <- lapply(ML_tree_rooted_clades, function(lst) lst[!gsub("_.*", "", lst) %in% uniprot_proteomes_euk_taxids])
prok_proteins_in_clades_raw_collapse <- unlist(lapply(prok_proteins_in_clades_raw, function(lst) paste0(lst, collapse=",")))
prok_groups_in_clades <- unlist(lapply(prok_proteins_in_clades_raw, function(lst) paste0(unique(gsub("_.*", "", lst)), collapse=",")))
n_prok_groups_in_clades <- 1 + str_count(prok_groups_in_clades, ",")
prok_proteins_in_clades <- lapply(prok_proteins_in_clades_raw, function(lst) sub("^[^_]*_", "", lst))
prok_species_in_clades <- unlist(lapply(prok_proteins_in_clades, function(lst) {
  paste0(sort(unique(unlist(strsplit(prokaryote_mmseqs2_clusters_taxids_in_tree$unique_member_seq_taxids[which(prokaryote_mmseqs2_clusters_taxids_in_tree$rep_seq %in% lst)], split=",")))), collapse=",")
}))
n_prok_species_in_clades <- unlist(lapply(prok_species_in_clades, function(lst) length(which(!is.na(unlist(strsplit(lst, split=",")))))))
n_bacteria_species_in_clades <- unlist(lapply(prok_species_in_clades, function(lst) length(which(unlist(strsplit(lst, split=",")) %in% bacteria_taxids))))
n_archaea_species_in_clades <- unlist(lapply(prok_species_in_clades, function(lst) length(which(unlist(strsplit(lst, split=",")) %in% archaea_taxids))))

n_species_in_clades <- n_euk_species_in_clades + n_prok_species_in_clades

curr_df <- data.frame(n_species_in_clades, n_euk_species_in_clades, n_prok_species_in_clades, n_bacteria_species_in_clades, n_archaea_species_in_clades, n_prok_groups_in_clades, euk_species_in_clades, prok_species_in_clades, prok_groups_in_clades, prok_proteins_in_clades_raw_collapse, euk_proteins_in_clades_raw_collapse)
curr_df$frac_euk_in_clades <- curr_df$n_euk_species_in_clades / curr_df$n_species_in_clades
curr_df$frac_prok_in_clades <- curr_df$n_prok_species_in_clades / curr_df$n_species_in_clades

n_euk_species_largest_euk_clade <- max(curr_df$n_euk_species_in_clades[curr_df$frac_euk_in_clades >= clade_purity_threshold])

if (sum(curr_df$frac_prok_in_clades >= clade_purity_threshold) > 0) {
  # Have pure prokaryotic clades
  curr_df$n_large_prokgroups_in_clade <- 0
  row_index_pure <- which(curr_df$frac_prok_in_clades >= clade_purity_threshold)
  for (i in row_index_pure) {
    curr_df$n_large_prokgroups_in_clade[i] <- count_large_prokgroups_in_clade(curr_df$prok_species_in_clades[i])
  }
  n_large_prokgroups_largest_prok_clade <- max(curr_df$n_large_prokgroups_in_clade)
  curr_df_pure_largest_prok_clade <- curr_df[which(curr_df$n_large_prokgroups_in_clade == n_large_prokgroups_largest_prok_clade),]
  
  n_prok_species_largest_prok_clade <- max(curr_df_pure_largest_prok_clade$n_prok_species_in_clades)
  n_prok_groups_largest_prok_clade <- max(n_prok_groups_in_clades)
  
  # Find largest prok clade by number of species prefiltered by largest number of prok groups
  prok_species_largest_prok_clade <- paste0(unique(unlist(strsplit(curr_df_pure_largest_prok_clade$prok_species_in_clades[which(curr_df_pure_largest_prok_clade$n_prok_species_in_clades == n_prok_species_largest_prok_clade)], split=","))), collapse=",")
  prok_groups_largest_prok_clade <- paste0(unique(unlist(strsplit(curr_df_pure_largest_prok_clade$prok_groups_in_clades[which(curr_df_pure_largest_prok_clade$n_prok_species_in_clades == n_prok_species_largest_prok_clade)], split=","))), collapse=",")
  prok_proteins_largest_prok_clade <- paste0(unique(unlist(strsplit(curr_df_pure_largest_prok_clade$prok_proteins_in_clades_raw_collapse[which(curr_df_pure_largest_prok_clade$n_prok_species_in_clades == n_prok_species_largest_prok_clade)], split=","))), collapse=",")
  
} else {
  # No pure prokaryotic clades
  n_large_prokgroups_largest_prok_clade <- 0
  n_prok_species_largest_prok_clade <- 0
  n_prok_groups_largest_prok_clade <- 0
  prok_species_largest_prok_clade <- ""
  prok_groups_largest_prok_clade <- ""
  prok_proteins_largest_prok_clade <- ""
}

euk_species_largest_euk_clade <- paste0(unique(unlist(strsplit(curr_df$euk_species_in_clades[which(curr_df$n_euk_species_in_clades == n_euk_species_largest_euk_clade)], split=","))), collapse=",")
euk_proteins_largest_euk_clade <- unique(unlist(strsplit(curr_df$euk_proteins_in_clades_raw_collapse[which(curr_df$n_euk_species_in_clades == n_euk_species_largest_euk_clade)], split=",")))

largest_monophyletic_clades <- data.frame(n_euk_species_largest_euk_clade, n_prok_species_largest_prok_clade, euk_species_largest_euk_clade, prok_species_largest_prok_clade, prok_groups_largest_prok_clade, prok_proteins_largest_prok_clade, n_large_prokgroups_largest_prok_clade)


prok_monophyletic_taxids <- unlist(strsplit(largest_monophyletic_clades$prok_groups_largest_prok_clade, split=","))
euk_monophyletic_taxids <- unlist(strsplit(largest_monophyletic_clades$euk_species_largest_euk_clade, split=","))

all_prok_proteins <- ML_tree$tip.label[!gsub("_.*", "", ML_tree$tip.label) %in% uniprot_proteomes_euk_taxids]

largest_monophyletic_clades$n_bacteria_species_largest_prok_clade <- sum(unique(unlist(strsplit(largest_monophyletic_clades$prok_species_largest_prok_clade, split=","))) %in% bacteria_taxids)
largest_monophyletic_clades$n_archaea_species_largest_prok_clade <- sum(unique(unlist(strsplit(largest_monophyletic_clades$prok_species_largest_prok_clade, split=","))) %in% archaea_taxids)

# Total species count
prok_species_df <- data.frame(prok_species_in_clades)
prok_species_df_seprows <- prok_species_df %>% separate_rows(prok_species_in_clades, sep=",") %>% filter(prok_species_in_clades != "")
n_total_prok_species <- length(unique(prok_species_df_seprows$prok_species_in_clades))
species_list <- unique(gsub("_.*", "", ML_tree$tip.label))
n_total_euk_species <- sum(species_list %in% uniprot_proteomes_euk_taxids)
frac_euk_species <- n_total_euk_species / (n_total_euk_species + n_total_prok_species)

# Determine origin
if (any(largest_monophyletic_clades$n_large_prokgroups_largest_prok_clade > 0 & largest_monophyletic_clades$n_prok_species_largest_prok_clade >= min_prok_species_in_clade_to_not_prune)) {
  origin_domain <- "Prokaryote"
} else if (frac_euk_species > min_frac_euk_species_for_euk_origin) {
  origin_domain <- "Eukaryote"
} else {
  origin_domain <- "Indeterminate"
}

if (bool_verbose) {
  print(paste0("Inferred origin: ", origin_domain))
  write.table(curr_df, paste0(out_dir, "/", MSA_id, "_clades.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
}

if (origin_domain == "Prokaryote") {
  ## Prokaryotic origin
  LCA_index <- get_mrca_of_set(species_tree, prok_monophyletic_taxids)
  
  # Drop small prokaryote clade tips from tree.
  curr_df_pure <- curr_df[curr_df$frac_prok_in_clades >= clade_purity_threshold,]
  prok_proteins_in_large_clades <- unique(unlist(strsplit(curr_df_pure$prok_proteins_in_clades_raw_collapse[which(curr_df_pure$n_large_prokgroups_in_clade > 0 | curr_df_pure$n_prok_species_in_clades >= min_prok_species_in_clade_to_not_prune)], split=",")))
  prok_proteins_in_small_clades <- all_prok_proteins[!all_prok_proteins %in% prok_proteins_in_large_clades]
  prok_proteins_to_drop <- prok_proteins_in_small_clades
  if (bool_verbose) {
    print(paste0("Dropping small prokaryote clade tips: ", paste0(prok_proteins_to_drop, collapse=",")))
  }
  ML_tree <- drop.tip(ML_tree, prok_proteins_to_drop)
  for (i in 1:length(ufboot_trees)) {
    ufboot_trees[[i]] <- drop.tip(ufboot_trees[[i]], prok_proteins_to_drop)
  }
  
  prok_proteins_largest_prok_clade_list <- unlist(strsplit(largest_monophyletic_clades$prok_proteins_largest_prok_clade, split=","))
  
  extract_subtree_for_clade <- function(tree, clade_list) {
    # Extract largest prok clade with only prokaryotes
    tree_largest_prok_clade <- get_subtree_with_tips(tree, only_tips=clade_list)
    tree_largest_prok_clade_subtree <- tree_largest_prok_clade$subtree
    
    if (is.rooted(tree_largest_prok_clade_subtree) & Ntip(tree_largest_prok_clade_subtree) > 2) {
      tree_largest_prok_clade_subtree <- unroot(tree_largest_prok_clade_subtree)
    }
    tree_largest_prok_clade_subtree$node.label <- NULL
    return(tree_largest_prok_clade_subtree)
  }
  
  # MAD rooting within the largest prokaryotic clade
  root_tree_from_rooted_subtree <- function(tree, tree_largest_prok_clade_subtree_rooted) {
    # Get the two monophyletic clades descending from the root of subtree
    root_daughter_nodes <- tree_largest_prok_clade_subtree_rooted$edge[which(tree_largest_prok_clade_subtree_rooted$edge[,1] == find_root(tree_largest_prok_clade_subtree_rooted)), 2]
    if (root_daughter_nodes[1] > Ntip(tree_largest_prok_clade_subtree_rooted)) {
      tree_largest_prok_clade_subtree_rooted_clade_subtree <- get_subtree_at_node(tree_largest_prok_clade_subtree_rooted, root_daughter_nodes[1] - Ntip(tree_largest_prok_clade_subtree_rooted))$subtree
      root_monophyletic_clade_tiplabels_daughter1 <- tree_largest_prok_clade_subtree_rooted_clade_subtree$tip.label
    } else {
      root_monophyletic_clade_tiplabels_daughter1 <- tree_largest_prok_clade_subtree_rooted$tip.label[root_daughter_nodes[1]]
    }
    if (root_daughter_nodes[2] > Ntip(tree_largest_prok_clade_subtree_rooted)) {
      tree_largest_prok_clade_subtree_rooted_clade_subtree <- get_subtree_at_node(tree_largest_prok_clade_subtree_rooted, root_daughter_nodes[2] - Ntip(tree_largest_prok_clade_subtree_rooted))$subtree
      root_monophyletic_clade_tiplabels_daughter2 <- tree_largest_prok_clade_subtree_rooted_clade_subtree$tip.label
    } else {
      root_monophyletic_clade_tiplabels_daughter2 <- tree_largest_prok_clade_subtree_rooted$tip.label[root_daughter_nodes[2]]
    }
    # Take the largest clade
    if (length(root_monophyletic_clade_tiplabels_daughter1) >= length(root_monophyletic_clade_tiplabels_daughter2)) {
      root_monophyletic_clade_tiplabels <- root_monophyletic_clade_tiplabels_daughter1
    } else {
      root_monophyletic_clade_tiplabels <- root_monophyletic_clade_tiplabels_daughter2
    }
    
    # Filter to find the root daughter clade that has only prokaryotes, important if largest prok clade is not monophyletic
    n_euk_species_in_clades <- unlist(lapply(ML_tree_rooted_clades_species, function(lst) sum(lst %in% uniprot_proteomes_euk_taxids)))
    n_prok_root_species_in_clades <- unlist(lapply(ML_tree_rooted_clades, function(lst) sum(lst %in% root_monophyletic_clade_tiplabels)))
    curr_clades_df <- data.frame(clade_index = 1:length(ML_tree_rooted_clades), n_euk_species_in_clades, n_prok_root_species_in_clades)
    curr_clades_df <- curr_clades_df[order(1/curr_clades_df$n_euk_species_in_clades, curr_clades_df$n_prok_root_species_in_clades, decreasing=TRUE),]
    root_monophyletic_clade_tiplabels <- ML_tree_rooted_clades[[curr_clades_df$clade_index[1]]]
    
    # Diagnostics, write out rooted original subtree
    if (bool_verbose) {
      print(paste0("Root clade: ", paste0(root_monophyletic_clade_tiplabels, collapse=", ")))
    }
    
    # Map subtree tips to the old tree for rooting
    original_root_index <- get_mrca_of_set(tree, root_monophyletic_clade_tiplabels)
    
    # If root is at the original default root, cannot directly root due to ambiguity when resolving a node in an unrooted tree with no explicit outgroup
    if (original_root_index == find_root(tree)) {
      if (bool_verbose) {
        print("New root same as default root. Rerooting...")
      }
      # First reroot at a eukaryotic tip, then reroot in prokaryotes
      tree <- root(tree, outgroup = euk_proteins_largest_euk_clade[1], resolve.root = TRUE)
      original_root_index <- get_mrca_of_set(tree, root_monophyletic_clade_tiplabels)
    }
    
    # Root the tree
    if (original_root_index > Ntip(tree)) {
      # Root is at a node
      tree <- root(tree, node = original_root_index, resolve.root = TRUE)
    } else {
      # Root is at a tip
      tree <- root(tree, outgroup = original_root_index, resolve.root = TRUE)
    }
    
    return(tree)
  }
  
  if (length(prok_proteins_largest_prok_clade_list) <= 2) {
    # Largest prokaryote clade has just 1 or 2 tips. Cannot root using MAD
    if (bool_verbose) {
      print(paste0("Outgroup rooting ML tree (largest prokaryote clade has 1-2 tips)"))
    }
    
    # Find root that achieves monophyly for largest prok clade
    root_index <- find_root_of_monophyletic_tips(ML_tree, prok_proteins_largest_prok_clade_list, as_MRCA=TRUE, is_rooted=FALSE)
    if (is.na(root_index)) {
      # If there's no root that achieves monophyly for largest prok clade, then use the LCA of largest prok clade
      # Reroot at a eukaryote monophyletic tip, then reroot
      tree <- root(ML_tree, outgroup = euk_proteins_largest_euk_clade[1], resolve.root = TRUE)
      root_index <- get_mrca_of_set(ML_tree, prok_proteins_largest_prok_clade_list)
    }
    
    if (root_index == find_root(ML_tree)) {
      # If root is at the original default root, cannot directly root due to ambiguity when resolving a node in an unrooted tree with no explicit outgroup
      # Root at randomly sampled tip in largest prok clade
      prok_proteins_largest_prok_clade_sampled_tip <- prok_proteins_largest_prok_clade_list[sample.int(length(prok_proteins_largest_prok_clade_list), size=1)]
      if (bool_verbose) {
        print(paste0("Rooted at randomly sampled tip in largest prok clade: ", prok_proteins_largest_prok_clade_sampled_tip))
      }
      ML_tree_rooted <- root(ML_tree, outgroup = prok_proteins_largest_prok_clade_sampled_tip, resolve.root = TRUE)
    } else if (root_index > Ntip(ML_tree)) {
      ML_tree_rooted <- root(ML_tree, node = root_index, resolve.root=TRUE)
    } else {
      ML_tree_rooted <- root(ML_tree, outgroup = root_index, resolve.root=TRUE)
    }
    write.tree(ML_tree_rooted, paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix, ".rooted"))
    
    if (bool_verbose) {
      print(paste0("Outgroup rooting UFBOOT trees (largest prokaryote clade has 1-2 tips)"))
    }
    ufboot_trees_rooted <- ufboot_trees # initialize
    for (i in 1:length(ufboot_trees)) {
      # Find root that achieves monophyly for largest prok clade
      root_index <- find_root_of_monophyletic_tips(ufboot_trees[[i]], prok_proteins_largest_prok_clade_list, as_MRCA=TRUE, is_rooted=FALSE)
      if (is.na(root_index)) {
        # If there's no root that achieves monophyly for largest prok clade, then use the LCA of largest prok clade
        # Reroot at a eukaryote monophyletic tip, then reroot
        tree <- root(ufboot_trees[[i]], outgroup = euk_proteins_largest_euk_clade[1], resolve.root = TRUE)
        root_index <- get_mrca_of_set(ufboot_trees[[i]], prok_proteins_largest_prok_clade_list)
      }
      
      if (root_index == find_root(ufboot_trees[[i]])) {
        # If root is at the original default root, cannot directly root due to ambiguity when resolving a node in an unrooted tree with no explicit outgroup
        # Root at randomly sampled tip in largest prok clade
        prok_proteins_largest_prok_clade_sampled_tip <- prok_proteins_largest_prok_clade_list[sample.int(length(prok_proteins_largest_prok_clade_list), size=1)]
        if (bool_verbose) {
          print(paste0("Rooted at randomly sampled tip in largest prok clade: ", prok_proteins_largest_prok_clade_sampled_tip))
        }
        ufboot_trees_rooted[[i]] <- root(ufboot_trees[[i]], outgroup = prok_proteins_largest_prok_clade_sampled_tip, resolve.root = TRUE)
      } else if (root_index > Ntip(ufboot_trees[[i]])) {
        ufboot_trees_rooted[[i]] <- root(ufboot_trees[[i]], node = root_index, resolve.root=TRUE)
      } else {
        ufboot_trees_rooted[[i]] <- root(ufboot_trees[[i]], outgroup = root_index, resolve.root=TRUE)
      }
    }
    write.tree(ufboot_trees_rooted, paste0(out_dir, "/", MSA_id, "_pruned.ufboot.rooted"))
  } else {
    # Largest prokaryote clade has more than 2 tips
    if (bool_verbose) {
      print(paste0("MAD rooting ML tree"))
    }
    tree_largest_prok_clade_subtree <- extract_subtree_for_clade(ML_tree, prok_proteins_largest_prok_clade_list)
    tree_largest_prok_clade_subtree_filename <- paste0(out_dir, "/", MSA_id, "_largest.prok.clade_pruned", ML_tree_suffix)
    write.tree(tree_largest_prok_clade_subtree, tree_largest_prok_clade_subtree_filename)
    outfile <- paste0(tree_largest_prok_clade_subtree_filename, ".rooted")
    if (!file.exists(outfile) | file.size(outfile) == 0) {
      system(paste0("mad -nt ", tree_largest_prok_clade_subtree_filename, " > /dev/null"))
    }
    tree_largest_prok_clade_subtree_rooted <- read.tree(paste0(tree_largest_prok_clade_subtree_filename, ".rooted"))
    ML_tree_rooted <- root_tree_from_rooted_subtree(ML_tree, tree_largest_prok_clade_subtree_rooted)
    write.tree(ML_tree_rooted, paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix, ".rooted"))
    
    if (bool_verbose) {
      print(paste0("MAD rooting UFBOOT trees"))
    }
    ufboot_trees_largest_prok_clade_subtrees <- ufboot_trees # initialize
    for (i in 1:length(ufboot_trees)) {
      ufboot_trees_largest_prok_clade_subtrees[[i]] <- extract_subtree_for_clade(ufboot_trees[[i]], prok_proteins_largest_prok_clade_list)
    }
    ufboot_trees_largest_prok_clade_subtree_filename <- paste0(out_dir, "/", MSA_id, "_largest.prok.clade_pruned.ufboot")
    write.tree(ufboot_trees_largest_prok_clade_subtrees, ufboot_trees_largest_prok_clade_subtree_filename)
    
    outfile <- paste0(ufboot_trees_largest_prok_clade_subtree_filename, ".rooted")
    if (file.exists(outfile)) {
      ufboot_trees_largest_prok_clade_subtree_rooted <- read.tree(outfile, keep.multi=TRUE)
      last_index_completed <- length(ufboot_trees_largest_prok_clade_subtree_rooted)
      if (last_index_completed < length(ufboot_trees_largest_prok_clade_subtrees)) {
        if (bool_verbose) {
          print(paste0("Resuming from prior run at ", last_index_completed, " of ", length(ufboot_trees_largest_prok_clade_subtrees)))
        }
        indexes_to_be_done <- (last_index_completed+1):length(ufboot_trees_largest_prok_clade_subtrees)
        ufboot_trees_largest_prok_clade_subtree_filename_tmp <- paste0(ufboot_trees_largest_prok_clade_subtree_filename, ".tmp")
        write.tree(ufboot_trees_largest_prok_clade_subtrees[indexes_to_be_done], ufboot_trees_largest_prok_clade_subtree_filename_tmp)
        
        system(paste0("mad -nt ", ufboot_trees_largest_prok_clade_subtree_filename_tmp, " > /dev/null"))
        
        system(paste0("cat ", outfile, " ", ufboot_trees_largest_prok_clade_subtree_filename_tmp, ".rooted > ", outfile, "_combined.tmp"))
        system(paste0("mv ", outfile, "_combined.tmp", " ", outfile))
        system(paste0("rm ", ufboot_trees_largest_prok_clade_subtree_filename_tmp))
        system(paste0("rm ", ufboot_trees_largest_prok_clade_subtree_filename_tmp, ".rooted"))
      }
    } else {
      system(paste0("mad -nt ", ufboot_trees_largest_prok_clade_subtree_filename, " > /dev/null"))
    }
    ufboot_trees_largest_prok_clade_subtree_rooted <- read.tree(paste0(ufboot_trees_largest_prok_clade_subtree_filename, ".rooted"))
    
    ufboot_trees_rooted <- ufboot_trees # initialize
    for (i in 1:length(ufboot_trees)) {
      ufboot_trees_rooted[[i]] <- root_tree_from_rooted_subtree(ufboot_trees[[i]], ufboot_trees_largest_prok_clade_subtree_rooted[[i]])
    }
    write.tree(ufboot_trees_rooted, paste0(out_dir, "/", MSA_id, "_pruned.ufboot.rooted"))
  }
  
} else if (origin_domain == "Indeterminate") {
  ## Indeterminate origin
  LCA_index <- get_mrca_of_set(species_tree, c(euk_monophyletic_taxids, prok_monophyletic_taxids))
  
  # Drop small prokaryote clade tips from tree.
  curr_df_pure <- curr_df[curr_df$frac_prok_in_clades >= clade_purity_threshold,]
  prok_proteins_in_large_clades <- unique(unlist(strsplit(curr_df_pure$prok_proteins_in_clades_raw_collapse[which(curr_df_pure$n_large_prokgroups_in_clade > 0 | curr_df_pure$n_prok_species_in_clades >= min_prok_species_in_clade_to_not_prune)], split=",")))
  prok_proteins_in_small_clades <- all_prok_proteins[!all_prok_proteins %in% prok_proteins_in_large_clades]
  prok_proteins_to_drop <- prok_proteins_in_small_clades
  if (bool_verbose) {
    print(paste0("Dropping small prokaryote clade tips: ", paste0(prok_proteins_to_drop, collapse=",")))
  }
  ML_tree <- drop.tip(ML_tree, prok_proteins_to_drop)
  for (i in 1:length(ufboot_trees)) {
    ufboot_trees[[i]] <- drop.tip(ufboot_trees[[i]], prok_proteins_to_drop)
  }
  
  # MAD rooting
  if (bool_verbose) {
    print(paste0("MAD rooting ML tree"))
  }
  write.tree(ML_tree, paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix))
  
  outfile <- paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix, ".rooted")
  if (!file.exists(outfile) | file.size(outfile) == 0) {
    system(paste0("mad -nt ", paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix), " > /dev/null"))
  }
  
  if (bool_verbose) {
    print(paste0("MAD rooting UFBOOT trees"))
  }
  
  ufboot_trees_filename <- paste0(out_dir, "/", MSA_id, "_pruned.ufboot")
  write.tree(ufboot_trees, ufboot_trees_filename)
  
  outfile <- paste0(ufboot_trees_filename, ".rooted")
  if (file.exists(outfile)) {
    ufboot_trees_rooted <- read.tree(outfile, keep.multi=TRUE)
    last_index_completed <- length(ufboot_trees_rooted)
    if (last_index_completed < length(ufboot_trees)) {
      if (bool_verbose) {
        print(paste0("Resuming from prior run at ", last_index_completed, " of ", length(ufboot_trees)))
      }
      indexes_to_be_done <- (last_index_completed+1):length(ufboot_trees)
      ufboot_trees_filename_tmp <- paste0(ufboot_trees_filename, ".tmp")
      write.tree(ufboot_trees[indexes_to_be_done], ufboot_trees_filename_tmp)
      system(paste0("mad -nt ", ufboot_trees_filename_tmp, " > /dev/null"))
      
      system(paste0("cat ", outfile, " ", ufboot_trees_filename_tmp, ".rooted > ", outfile, "_combined.tmp"))
      system(paste0("mv ", outfile, "_combined.tmp", " ", outfile))
      system(paste0("rm ", ufboot_trees_filename_tmp))
      system(paste0("rm ", ufboot_trees_filename_tmp, ".rooted"))
    }
    
  } else {
    system(paste0("mad -nt ", ufboot_trees_filename, " > /dev/null"))
  }
  
} else if (origin_domain == "Eukaryote") {
  ## Eukaryotic origin
  LCA_index <- get_mrca_of_set(species_tree, euk_monophyletic_taxids)
  
  # Drop all prokaryote clade tips from tree
  prok_proteins_to_drop <- all_prok_proteins
  if (bool_verbose) {
    print(paste0("Dropping all prokaryote clade tips: ", paste0(prok_proteins_to_drop, collapse=",")))
  }
  ML_tree <- drop.tip(ML_tree, prok_proteins_to_drop)
  for (i in 1:length(ufboot_trees)) {
    ufboot_trees[[i]] <- drop.tip(ufboot_trees[[i]], prok_proteins_to_drop)
  }
  
  # MAD rooting
  if (bool_verbose) {
    print(paste0("MAD rooting ML tree"))
  }
  write.tree(ML_tree, paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix))
  
  outfile <- paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix, ".rooted")
  if (!file.exists(outfile) | file.size(outfile) == 0) {
    system(paste0("mad -nt ", paste0(out_dir, "/", MSA_id, "_pruned", ML_tree_suffix), " > /dev/null"))
  }
  
  if (bool_verbose) {
    print(paste0("MAD rooting UFBOOT trees"))
  }
  
  ufboot_trees_filename <- paste0(out_dir, "/", MSA_id, "_pruned.ufboot")
  write.tree(ufboot_trees, ufboot_trees_filename)
  
  outfile <- paste0(ufboot_trees_filename, ".rooted")
  if (file.exists(outfile)) {
    ufboot_trees_rooted <- read.tree(outfile, keep.multi=TRUE)
    last_index_completed <- length(ufboot_trees_rooted)
    if (last_index_completed < length(ufboot_trees)) {
      if (bool_verbose) {
        print(paste0("Resuming from prior run at ", last_index_completed, " of ", length(ufboot_trees)))
      }
      indexes_to_be_done <- (last_index_completed+1):length(ufboot_trees)
      ufboot_trees_filename_tmp <- paste0(ufboot_trees_filename, ".tmp")
      write.tree(ufboot_trees[indexes_to_be_done], ufboot_trees_filename_tmp)
      system(paste0("mad -nt ", ufboot_trees_filename_tmp, " > /dev/null"))
      
      system(paste0("cat ", outfile, " ", ufboot_trees_filename_tmp, ".rooted > ", outfile, "_combined.tmp"))
      system(paste0("mv ", outfile, "_combined.tmp", " ", outfile))
      system(paste0("rm ", ufboot_trees_filename_tmp))
      system(paste0("rm ", ufboot_trees_filename_tmp, ".rooted"))
    }
    
  } else {
    system(paste0("mad -nt ", ufboot_trees_filename, " > /dev/null"))
  }
}

LCA_euks_index <- get_mrca_of_set(species_tree, euk_monophyletic_taxids)

# Get LCA node label
LCA_node_ALE <- species_tree_labels[LCA_index]
LCA_node_euks_ALE <- species_tree_labels[LCA_euks_index]
if (bool_verbose) {
  print(paste0("LCA node: ", LCA_node_ALE))
  print(paste0("LCA euk node: ", LCA_node_euks_ALE))
}

dropped_tips <- c(dup_proteins, remove_mtDNA_proteins, prok_proteins_to_drop)

LCA_df <- data.frame(MSA_id = MSA_id, n_euk_species_largest_euk_clade = largest_monophyletic_clades$n_euk_species_largest_euk_clade, n_prok_species_largest_prok_clade = largest_monophyletic_clades$n_prok_species_largest_prok_clade, LCA_node_ALE = LCA_node_ALE, LCA_node_euks_ALE = LCA_node_euks_ALE, origin_domain = origin_domain, dropped_tips = paste0(dropped_tips, collapse=","))

# Check for completion
ufboot_rooted_trees <- read.tree(paste0(out_dir, "/", MSA_id, "_pruned.ufboot.rooted"))
if (length(ufboot_rooted_trees) < n_ufboot_trees) {
  print(paste0(MSA_id, " task failed. Fewer than ", n_ufboot_trees, " ufboot trees successfully rooted."))
  quit(save="no")
}


## Write out
# write.table(LCA_df, here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

