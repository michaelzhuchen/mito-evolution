# Load libraries
library(tidyverse)
library(ape)
library(castor)
library(ggplot2)
library(here)

### Read in datasets
# Read in experimental and mtDNA mito proteins
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2025.09.30.tsv"), sep="\t", header=TRUE)

# Read in OG domain origins
origin_table <- read.table(here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t")
colnames(origin_table) <- c("OG_id", "n_euk_species_largest_euk_clade", "n_prok_species_largest_prok_clade", "LCA_node", "LCA_node_euks", "origin_domain", "dropped_tips")
origin_table$OG_id <- gsub(".faa_clipkit.gappy.msa", "", origin_table$OG_id, fixed=TRUE)
prok_origin_OG_ids <- origin_table$OG_id[origin_table$origin_domain == "Prokaryote"]

# Read in OGs
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")
ogs_long_primary <- ogs_long %>% filter(BOOL_PRIMARY_OG)

## Load and process DeepLoc results
# Read in DeepLoc results
deeploc_results <- read.table(file.path(supplemental_data_directory, "TableS7_retrained_DeepLoc_predictions.tsv"), header=TRUE)
colnames(deeploc_results) <- c("Protein_ID", "Mitochondrion")
deeploc_results$taxid <- gsub("_.*", "", deeploc_results$Protein_ID)
# Fix organelle-encoded proteins
all_nonmito_accessions <- read.table(here("data/deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
all_mito_accessions <- read.table(here("data/deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1
deeploc_results$Mitochondrion[deeploc_results$Protein_ID %in% all_nonmito_accessions] <- 0
deeploc_results$Mitochondrion[deeploc_results$Protein_ID %in% all_mito_accessions] <- 1
# Mark experimentally-defined proteomes
completed_mitoproteomes_species_list <- c("9606", "559292", "3702", "1257118", "5689", "185431", "5741", "32595")
deeploc_results$Mitochondrion[deeploc_results$taxid %in% completed_mitoproteomes_species_list & !deeploc_results$Protein_ID %in% gold_gene_accession_OG_id_df$gene_accession] <- 0
deeploc_results$Mitochondrion[deeploc_results$Protein_ID %in% gold_gene_accession_OG_id_df$gene_accession] <- 1

# Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

eukaryote_reference_species_list <- completed_mitoproteomes_species_list

# Read in species tree and select key ancestral nodes connecting focal species
species_tree <- read.tree(here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample_v2.contree"))
tree <- keep.tip(species_tree, eukaryote_reference_species_list)

## Set mito localization thresholds for Mk model and parsimony model
deeploc_thresholds <- read.csv(here("data/deeploc", "deeploc_thresholds.csv"))
deeploc_mito_threshold <- deeploc_thresholds$threshold[deeploc_thresholds$label == "Mitochondrion"]
mito_localization_prob_mk_threshold <- deeploc_mito_threshold
mito_localization_prob_parsimony_threshold <- 0.5

# Filter DeepLocMC predictions by minimum number of species with predicted mito protein
n_mito_species_threshold <- 5
deeploc_results_ogs <- merge(deeploc_results, ogs_long, by.x=c("Protein_ID", "taxid"), by.y=c("accession", "taxid"))
deeploc_results_exclude <- deeploc_results_ogs %>% filter(Mitochondrion >= deeploc_mito_threshold) %>% group_by(Orthogroup) %>% filter(length(unique(taxid)) < n_mito_species_threshold)
deeploc_results_exclude_OG_ids <- unique(deeploc_results_exclude$Orthogroup)
deeploc_results_exclude_OG_ids <- deeploc_results_exclude_OG_ids[!deeploc_results_exclude_OG_ids %in% gold_gene_accession_OG_id_df$OG_id]
## Write out OGs that have few species with predicted mito proteins
# write.table(deeploc_results_exclude_OG_ids, here('data/reconstruction', 'deeploc_results_exclude_lt5mitospp_OG_ids.txt'), sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)

tax_levels_df <- data.frame(label = tree$node.label)
tree_labels <- c(tree$tip.label, tree$node.label)
tax_levels_df$node_index <- match(tax_levels_df$label, tree_labels)

## Read in Eukaryota parent PhROGs
n_minimum_species_leca <- 4
parent_progs_long_raw <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(parent_progs_long_raw) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
parent_progs_long <- parent_progs_long_raw %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))
parent_progs_long_primary <- parent_progs_long_raw %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% filter(BOOL_primary_OG) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))

# Filter by mito localization prob by Mk or parsimony
parent_progs_long_primary_mito <- parent_progs_long_primary %>% group_by(PROG_id) %>% filter(mito_localization_prob_mk >= mito_localization_prob_mk_threshold | mito_localization_prob_parsimony >= mito_localization_prob_parsimony_threshold)
# Filter out DeepLocMC predictions below a minimum number of species with predicted mito protein
parent_progs_long_primary_mito <- parent_progs_long_primary_mito %>% filter(!OG_id %in% deeploc_results_exclude_OG_ids)
parent_mito_PhROG_ids <- unique(parent_progs_long_primary_mito$PROG_id)
## Write out parent mito PhROG ids
# write.table(parent_mito_PhROG_ids, here("data/reconstruction", "parent_mito_PhROG_ids_corespecies_deeplocmc.min5spp_Eukaryota.parent.2supergroups.min4spp.filter.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# Get Eukaryota_parent
parent_progs_long_primary_mito_eukaryota <- parent_progs_long_primary_mito %>% filter(label == "Node34_Eukaryota")
parent_progs_long_primary_mito_eukaryota <- parent_progs_long_primary_mito_eukaryota %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)
eukaryota_parent_mito_PhROG_ids <- unique(parent_progs_long_primary_mito_eukaryota$PROG_id)
parent_progs_long_eukaryota <- parent_progs_long %>% filter(label == "Node34_Eukaryota")

# Filter for prok origin for HGT detection
parent_progs_long_prokorigin <- parent_progs_long %>% filter(OG_id %in% prok_origin_OG_ids)

## Get PhROGs at each taxonomic level, traversing from root to leaves
progs_long_agg <- c()
for (i in 1:nrow(tax_levels_df)) {
  tax_levels_df_curr <- tax_levels_df[i,]
  print(tax_levels_df_curr$label)
  
  progs_long <- read.table(file.path("data/phylogenetically_resolved_orthogroups", "PhROGs_long", paste0("PhROGs_at_", tax_levels_df_curr$label, "_long.tsv")), sep="\t", header=TRUE)
  colnames(progs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
  
  progs_long <- progs_long %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id))
  
  # For NA mito localization probs, retrieve the mito localization prob from the most recent parent
  if (tax_levels_df_curr$label != tree$node.label[1]) {
    # If not the root, impute from the most recent parent
    ancestor_node_index <- tree$edge[tree$edge[,2] == tax_levels_df_curr$node_index,1]
    ancestor_node_label <- tree_labels[ancestor_node_index]
    progs_long_agg_parent_curr <- progs_long_agg %>% filter(grepl(ancestor_node_label, PROG_id, fixed=TRUE))
    parent_OG_protein_ids <- paste0(progs_long_agg_parent_curr$OG_id, "_", progs_long_agg_parent_curr$protein_id)
    
    # Match by OG_proteinID to ensure correct matching of primary/nonprimary OGs
    missing_index <- which(is.na(progs_long$mito_localization_prob_mk))
    missing_child_OG_protein_ids <- paste0(progs_long$OG_id[missing_index], "_", progs_long$protein_id[missing_index])
    progs_long$mito_localization_prob_mk[missing_index] <- progs_long_agg_parent_curr$mito_localization_prob_mk[match(missing_child_OG_protein_ids, parent_OG_protein_ids)]
    
    missing_index <- which(is.na(progs_long$mito_localization_prob_parsimony))
    missing_child_OG_protein_ids <- paste0(progs_long$OG_id[missing_index], "_", progs_long$protein_id[missing_index])
    progs_long$mito_localization_prob_parsimony[missing_index] <- progs_long_agg_parent_curr$mito_localization_prob_parsimony[match(missing_child_OG_protein_ids, parent_OG_protein_ids)]
    
  } else {
    # If the root, impute from Eukaryota_parent
    parent_OG_protein_ids <- paste0(parent_progs_long$OG_id, "_", parent_progs_long$protein_id)
    
    # Match by OG_proteinID to ensure correct matching of primary/nonprimary OGs
    missing_index <- which(is.na(progs_long$mito_localization_prob_mk))
    missing_child_OG_protein_ids <- paste0(progs_long$OG_id[missing_index], "_", progs_long$protein_id[missing_index])
    progs_long$mito_localization_prob_mk[missing_index] <- parent_progs_long$mito_localization_prob_mk[match(missing_child_OG_protein_ids, parent_OG_protein_ids)]
    
    missing_index <- which(is.na(progs_long$mito_localization_prob_parsimony))
    missing_child_OG_protein_ids <- paste0(progs_long$OG_id[missing_index], "_", progs_long$protein_id[missing_index])
    progs_long$mito_localization_prob_parsimony[missing_index] <- parent_progs_long$mito_localization_prob_parsimony[match(missing_child_OG_protein_ids, parent_OG_protein_ids)]
  }
  
  progs_long_agg <- rbind(progs_long_agg, progs_long)
  
}

### Generate gene table for mito homologs. Primary OGs only
gene_table <- c()
for (i in 1:nrow(tax_levels_df)) {
  tax_levels_df_curr <- tax_levels_df[i,]
  print(tax_levels_df_curr$label)
  
  progs_long_filter <- progs_long_agg %>% filter(grepl(tax_levels_df_curr$label, PROG_id, fixed=TRUE)) %>% filter(BOOL_primary_OG) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))
  
  # Filter by mito localization prob by Mk or parsimony
  progs_long_filter <- progs_long_filter %>% group_by(PROG_id) %>% filter(mito_localization_prob_mk >= mito_localization_prob_mk_threshold | mito_localization_prob_parsimony >= mito_localization_prob_parsimony_threshold)
  
  # Filter out DeepLocMC predictions below a minimum number of species with predicted mito protein
  progs_long_filter <- progs_long_filter %>% filter(!OG_id %in% deeploc_results_exclude_OG_ids)
  
  # Select PhROGs present at the tax level or inherited from earlier tax level
  if (tax_levels_df_curr$node_index != Ntip(tree)+1) {
    # If not the root, find rows that are present at the current tax level or were present in ancestor
    ancestor_node_index <- tree$edge[tree$edge[,2] == tax_levels_df_curr$node_index,1]
    gene_table_ancestor <- gene_table %>% filter(node == ancestor_node_index)
    
    parent_OG_protein_ids <- paste0(gene_table_ancestor$origin_family, "_", gene_table_ancestor$protein_id)
    OG_protein_ids <- paste0(progs_long_filter$OG_id, "_", progs_long_filter$protein_id)
    present_PhROG_ids <- unique(progs_long_filter$PROG_id[which(progs_long_filter$label == tax_levels_df_curr$label | OG_protein_ids %in% parent_OG_protein_ids)])
    progs_long_filter <- progs_long_filter %>% filter(PROG_id %in% present_PhROG_ids)
    
  } else {
    
    # If root, find rows that are present at the current tax level or in Eukaryota_parent mito
    parent_OG_protein_ids <- paste0(parent_progs_long_primary_mito_eukaryota$OG_id, "_", parent_progs_long_primary_mito_eukaryota$protein_id)
    OG_protein_ids <- paste0(progs_long_filter$OG_id, "_", progs_long_filter$protein_id)
    present_PhROG_ids <- unique(progs_long_filter$PROG_id[which(progs_long_filter$label == tax_levels_df_curr$label | OG_protein_ids %in% parent_OG_protein_ids)])
    progs_long_filter <- progs_long_filter %>% filter(PROG_id %in% present_PhROG_ids)
    
    # Require that it be present in a euk supergroup indicating breadth of coverage
    permitted_leca_duplication_LCA_nodelabels <- c("Eukaryota", "Opimoda", "Amorphea_CRuMs", "Amorphea", "Obazoa", "Diphoda", "Diaphorectickes", "CAM_Haptista", "CAM")
    progs_long_filter <- progs_long_filter %>% group_by(PROG_id) %>% filter(gsub("Node[0-9]+_", "", label) %in% permitted_leca_duplication_LCA_nodelabels)
    
    # Minimum species filter
    progs_long_filter <- progs_long_filter %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)
  }
  
  progs_long_filter <- progs_long_filter %>% rowwise() %>% mutate(mito_localization_prob_max = max(mito_localization_prob_mk, mito_localization_prob_parsimony, na.rm=TRUE))
  
  gene_table_curr <- data.frame(node = tax_levels_df_curr$node_index, gene_id = progs_long_filter$PROG_id, origin_family = progs_long_filter$OG_id, protein_id = progs_long_filter$protein_id, mito_localization_prob = progs_long_filter$mito_localization_prob_max, stringsAsFactors = FALSE)
  gene_table <- rbind(gene_table, gene_table_curr)
}

# Get species mito proteomes
for (curr_taxid in eukaryote_reference_species_list) {
  if (curr_taxid %in% completed_mitoproteomes_species_list) {
    # If MitoCarta species, use experimental data
    gold_gene_accession_OG_id_df_curr <- gold_gene_accession_OG_id_df %>% filter(taxid == curr_taxid)
    gene_table_curr <- data.frame(node = match(curr_taxid, tree_labels), gene_id = gold_gene_accession_OG_id_df_curr$gene_accession, origin_family = gold_gene_accession_OG_id_df_curr$OG_id, protein_id = gold_gene_accession_OG_id_df_curr$gene_accession, mito_localization_prob = 1)
  } else {
    # If not MitoCarta species, use deeploc data
    deeploc_results_curr_species <- deeploc_results %>% filter(Mitochondrion >= deeploc_mito_threshold) %>% filter(taxid == curr_taxid)
    deeploc_results_curr_species$OG_id <- ogs_long_primary$Orthogroup[match(deeploc_results_curr_species$Protein_ID, ogs_long_primary$accession)]
    deeploc_results_curr_species <- deeploc_results_curr_species %>% filter(!OG_id %in% deeploc_results_exclude_OG_ids)
    gene_table_curr <- data.frame(node = match(curr_taxid, tree_labels), gene_id = deeploc_results_curr_species$Protein_ID, origin_family = deeploc_results_curr_species$OG_id, protein_id = deeploc_results_curr_species$Protein_ID, mito_localization_prob = deeploc_results_curr_species$Mitochondrion)
  }
  
  gene_table <- rbind(gene_table, gene_table_curr)
}
rownames(gene_table) <- 1:nrow(gene_table)

### Add mtDNA proteins from mtDNA-rich eukaryotes as gold mito proteins
# Read in hmmsearch results from refined OGs to mtDNA proteins in mtDNA-rich eukaryotes
mtdna_rich_euks_hmmsearch_tophit <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_mtDNA_rich_euks_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"), header=FALSE)
colnames(mtdna_rich_euks_hmmsearch_tophit) <- c("gene_accession", "OG_id")
mtdna_rich_euks_hmmsearch_tophit$taxid <- gsub("_.*", "", mtdna_rich_euks_hmmsearch_tophit$gene_accession)
gold_gene_accession_OG_id_other_species_mtdna_mantamonas_df <- gold_gene_accession_OG_id_df %>% filter(grepl("^2983909_", gene_accession))
mtdna_rich_euks_proteins_df <- rbind(gold_gene_accession_OG_id_other_species_mtdna_mantamonas_df, mtdna_rich_euks_hmmsearch_tophit)
mtdna_rich_euks_proteins_df$origin <- "Prokaryote"
# Remove homology failure for rps6 and remove tRNA-W(CCA) gene that is not a protein
manual_remove_protein_ids <- c("2983909_2983909_complement_18442_18756__rps6", "2983909_2983909_31046_31118_trnW_cca_")
manual_remove_OG_ids <- gold_gene_accession_OG_id_df$OG_id[gold_gene_accession_OG_id_df$gene_accession %in% manual_remove_protein_ids]
mtdna_rich_euks_proteins_df <- mtdna_rich_euks_proteins_df %>% filter(!OG_id %in% manual_remove_OG_ids)

# Add mtdna rich euk proteins to the gene table
gene_table_add_mtdna_rich_euks_agg <- c()
for (curr_OG_id in unique(mtdna_rich_euks_proteins_df$OG_id)) {
  mtdna_rich_euks_proteins_df_curr <- mtdna_rich_euks_proteins_df %>% filter(OG_id == curr_OG_id)
  BOOL_present_jakobid <- any(mtdna_rich_euks_proteins_df_curr$taxid %in% c("48483", "505711"))
  BOOL_present_crums <- any(mtdna_rich_euks_proteins_df_curr$taxid %in% c("2983909"))
  eukaryota_node <- which(tree_labels == "Node34_Eukaryota")
  diphoda_node <- which(tree_labels == "Node35_Diphoda")
  opimoda_node <- which(tree_labels == "Node139_Opimoda")
  gene_accession_curr <- mtdna_rich_euks_proteins_df_curr$gene_accession
  gene_accession_curr_jakobid <- mtdna_rich_euks_proteins_df_curr$gene_accession[mtdna_rich_euks_proteins_df_curr$taxid %in% c("48483", "505711")]
  gene_accession_curr_crums <- mtdna_rich_euks_proteins_df_curr$gene_accession[mtdna_rich_euks_proteins_df_curr$taxid %in% c("2983909")]
  
  nodes_to_add <- c(eukaryota_node)
  if (BOOL_present_jakobid) {
    nodes_to_add <- c(nodes_to_add, diphoda_node)
  }
  if (BOOL_present_crums) {
    nodes_to_add <- c(nodes_to_add, opimoda_node)
  }
  
  # If it's already in the gene table, then fill in the intermediate nodes downstream of Diphoda/Opimoda
  if (curr_OG_id %in% gene_table$origin_family) {
    gene_table_curr <- gene_table %>% filter(origin_family == curr_OG_id)
    gene_accession_curr <- c(gene_accession_curr, gene_table_curr$protein_id)
    present_nodes <- unique(gene_table_curr$node)
    present_labels <- tree_labels[present_nodes]
    
    if (BOOL_present_jakobid) {
      subtree <- get_subtree_at_node(tree, diphoda_node - Ntip(tree))$subtree
      present_nodes_diphoda <- present_nodes[present_labels %in% c(subtree$tip.label, subtree$node.label)]
      for (present_nodes_diphoda_curr in present_nodes_diphoda) {
        curr_node <- present_nodes_diphoda_curr
        n_splits <- 1
        intermediate_nodes <- c()
        while (curr_node != diphoda_node) {
          curr_node <- get_ancestral_nodes(tree, curr_node, Nsplits=1) + Ntip(tree)
          intermediate_nodes <- c(intermediate_nodes, curr_node)
        }
        nodes_to_add <- c(nodes_to_add, intermediate_nodes)
      }
    }
    if (BOOL_present_crums) {
      subtree <- get_subtree_at_node(tree, opimoda_node - Ntip(tree))$subtree
      present_nodes_opimoda <- present_nodes[present_labels %in% c(subtree$tip.label, subtree$node.label)]
      for (present_nodes_opimoda_curr in present_nodes_opimoda) {
        curr_node <- present_nodes_opimoda_curr
        n_splits <- 1
        intermediate_nodes <- c()
        while (curr_node != opimoda_node) {
          curr_node <- get_ancestral_nodes(tree, curr_node, Nsplits=1) + Ntip(tree)
          intermediate_nodes <- c(intermediate_nodes, curr_node)
        }
        nodes_to_add <- c(nodes_to_add, intermediate_nodes)
      }
    }
    nodes_to_add <- unique(nodes_to_add)
  }
  
  # Add new entries to gene_table
  gene_table_add <- c()
  for (curr_node in nodes_to_add) {
    gene_accession_curr_node <- c()
    gene_table_curr_index <- which(gene_table$node == curr_node & gene_table$origin_family == curr_OG_id)
    # Add if not found in the existing gene table
    if (length(gene_table_curr_index) == 0) {
      full_tree_node <- which(species_tree$node.label == tree_labels[curr_node])
      subtree <- get_subtree_at_node(species_tree, full_tree_node)$subtree
      gene_accession_curr_node <- gene_accession_curr[which(gsub("_.*", "", gene_accession_curr) %in% subtree$tip.label)]
      if ("Node35_Diphoda" %in% subtree$node.label) {
        gene_accession_curr_node <- c(gene_accession_curr_node, gene_accession_curr_jakobid)
      }
      if ("Node139_Opimoda" %in% subtree$node.label) {
        gene_accession_curr_node <- c(gene_accession_curr_node, gene_accession_curr_crums)
      }
      
      gene_table_add_curr <- data.frame(node = curr_node, gene_id = paste0(curr_OG_id, "_", tree_labels[curr_node], "_PhROG0000001"), origin_family = curr_OG_id, protein_id = gene_accession_curr_node, mito_localization_prob = 1)
      gene_table_add <- rbind(gene_table_add, gene_table_add_curr)
    }
  }
  
  gene_table_add_mtdna_rich_euks_agg <- rbind(gene_table_add_mtdna_rich_euks_agg, gene_table_add)
}
gene_table <- rbind(gene_table, gene_table_add_mtdna_rich_euks_agg)
gene_table <- gene_table[!duplicated(gene_table),]

### Identify vertically inherited, gained, and lost proteins on each branch
edge_list <- data.frame(tree$edge)
colnames(edge_list) <- c("parent", "child")

# Map node numbers to labels for tips
tip_labels <- tree$tip.label
node_labels <- as.character((length(tip_labels) + 1):(length(tip_labels) + tree$Nnode))
all_labels <- c(tip_labels, node_labels)
names(all_labels) <- as.character(1:(length(tip_labels) + tree$Nnode))

edge_flows <- list()
edge_flows_gains <- list()
edge_flows_losses <- list()
edge_flows_nonvertical <- list()

multiple_mappings <- c()
inexact_matches_agg <- c()

for (i in 1:nrow(edge_list)) {
  parent <- edge_list[i, "parent"]
  child <- edge_list[i, "child"]
  
  parent_id <- parent
  child_id <- child
  
  genes_p <- gene_table %>% filter(node == parent_id)
  genes_c <- gene_table %>% filter(node == child_id)
  
  # Join on exact proteins
  shared_fam_temp <- inner_join(genes_p, genes_c, by = c("origin_family", "protein_id"))
  colnames(shared_fam_temp)[which(colnames(shared_fam_temp) == "mito_localization_prob.y")] <- "mito_localization_prob"
  
  # Find cases where a single daughter PhROG maps to multiple parent PhROGs causing inconsistent counts
  multiple_mapping <- shared_fam_temp %>% group_by(origin_family, gene_id.y) %>% summarize(node = child_id, n_parent_PhROGs = length(unique(gene_id.x)), parent_PhROGs = paste0(sort(unique(gene_id.x)), collapse=","), proteins = paste0(sort(unique(protein_id)), collapse=","), .groups = "drop") %>% filter(n_parent_PhROGs > 1)
  multiple_mappings <- rbind(multiple_mappings, multiple_mapping)
  
  # Find cases where a single daughter PhROG has some members mapping to multiple parent PhROGs causing inconsistent counts
  genes_c_absent_p <- genes_c[!genes_c$protein_id %in% genes_p$protein_id,]
  shared_fam_temp_inexact_match <- genes_c_absent_p[genes_c_absent_p$gene_id %in% shared_fam_temp$gene_id.y,]
  inexact_matches <- shared_fam_temp[match(shared_fam_temp_inexact_match$gene_id, shared_fam_temp$gene_id.y),]
  inexact_matches$protein_id <- shared_fam_temp_inexact_match$protein_id
  inexact_matches_agg <- rbind(inexact_matches_agg, inexact_matches)
  shared_fam_temp <- rbind(shared_fam_temp, inexact_matches)
  
  shared_fam <- shared_fam_temp %>% group_by(gene_id.y) %>% summarize(origin_family = unique(origin_family), protein_id = paste0(sort(unique(protein_id)), collapse=","), mito_localization_prob = mean(mito_localization_prob)) # keep just one row per child PhROG
  
  # Identify gains
  gains_fam <- genes_c[!genes_c$gene_id %in% shared_fam$gene_id.y,]
  gains_fam <- gains_fam %>% group_by(gene_id) %>% summarize(origin_family = unique(origin_family), protein_id = paste0(sort(unique(protein_id)), collapse=","), mito_localization_prob = mean(mito_localization_prob))
  
  # Identify losses
  losses_fam <- genes_p[!genes_p$gene_id %in% genes_p$gene_id[genes_p$protein_id %in% genes_c$protein_id],]
  
  colnames(shared_fam_temp)[1] <- "node"
  colnames(shared_fam_temp)[2] <- "gene_id"
  losses_fam <- rbind(losses_fam, shared_fam_temp[,c("node", "gene_id", "origin_family", "protein_id", "mito_localization_prob")])
  
  losses_fam <- losses_fam %>% group_by(gene_id) %>% summarize(origin_family = unique(origin_family), protein_id = paste0(sort(unique(protein_id)), collapse=","), mito_localization_prob = mean(mito_localization_prob))
  
  # Find non-parsimonious discordances: OGs that have PhROGs lost and regained on the same branch.
  overlap_fam <- intersect(gains_fam$origin_family, losses_fam$origin_family)
  overlap_fam <- gains_fam[gains_fam$origin_family %in% overlap_fam,]
  
  flow_summary <- shared_fam %>%
    group_by(origin_family, protein_id) %>% summarize(count = n(), mito_localization_prob = mean(mito_localization_prob), .groups = "keep") %>%
    mutate(parent = child_id+0.75, child = child_id)
  
  flow_gains_summary <- gains_fam %>%
    group_by(origin_family, protein_id) %>% summarize(count = n(), mito_localization_prob = mean(mito_localization_prob), .groups = "keep") %>%
    mutate(parent = child_id+0.25, child = child_id)
  
  flow_losses_summary <- losses_fam %>%
    group_by(origin_family, protein_id) %>% summarize(count = n(), mito_localization_prob = mean(mito_localization_prob), .groups = "keep") %>%
    mutate(parent = parent_id, child = child_id+0.75)
  
  flow_nonvertical_summary <- overlap_fam %>%
    group_by(origin_family, protein_id) %>% summarize(count = n(), mito_localization_prob = mean(mito_localization_prob), .groups = "keep") %>%
    mutate(parent = child_id+0.75, child = child_id+0.25)
  
  edge_flows[[i]] <- flow_summary
  edge_flows_gains[[i]] <- flow_gains_summary
  edge_flows_losses[[i]] <- flow_losses_summary
  edge_flows_nonvertical[[i]] <- flow_nonvertical_summary
}

# Add the root OGs with prokaryotic origin
root_id <- find_root(tree)
root_fam <- gene_table %>% filter(node == root_id) %>% filter(origin_family %in% prok_origin_OG_ids)
root_fam <- root_fam %>% group_by(gene_id) %>% summarize(origin_family = unique(origin_family), protein_id = paste0(sort(unique(protein_id)), collapse=","), mito_localization_prob = mean(mito_localization_prob))
flow_root_prok_summary <- root_fam %>%
  group_by(origin_family, protein_id) %>% summarize(count = n(), mito_localization_prob = mean(mito_localization_prob), .groups = "keep") %>%
  mutate(parent = root_id+0.25, child = root_id)
# Add the root OGs with non-prokaryotic origin
root_id <- find_root(tree)
root_fam <- gene_table %>% filter(node == root_id) %>% filter(!origin_family %in% prok_origin_OG_ids)
root_fam <- root_fam %>% group_by(gene_id) %>% summarize(origin_family = unique(origin_family), protein_id = paste0(sort(unique(protein_id)), collapse=","), mito_localization_prob = mean(mito_localization_prob))
flow_root_nonprok_summary <- root_fam %>%
  group_by(origin_family, protein_id) %>% summarize(count = n(), mito_localization_prob = mean(mito_localization_prob), .groups = "keep") %>%
  mutate(parent = root_id+0.75, child = root_id)

# Combine
edge_vertical_df <- bind_rows(edge_flows)
edge_gains_df <- bind_rows(edge_flows_gains)
edge_losses_df <- bind_rows(edge_flows_losses)
edge_nonvertical_df <- bind_rows(edge_flows_nonvertical)
edge_root_prok_df <- bind_rows(flow_root_prok_summary)
edge_root_nonprok_df <- bind_rows(flow_root_nonprok_summary)
edge_vertical_df$is_singleton <- 0
edge_gains_df$is_singleton <- 1
edge_losses_df$is_singleton <- 0
edge_nonvertical_df$is_singleton <- 0
edge_root_prok_df$is_singleton <- 0
edge_root_nonprok_df$is_singleton <- 0
edge_vertical_df$is_nonvertical <- 0
edge_gains_df$is_nonvertical <- 0
edge_losses_df$is_nonvertical <- 0
edge_nonvertical_df$is_nonvertical <- 1
edge_root_prok_df$is_nonvertical <- 0
edge_root_nonprok_df$is_nonvertical <- 0
edge_root_df <- rbind(edge_root_prok_df, edge_root_nonprok_df)
edge_df <- rbind(edge_vertical_df, edge_gains_df, edge_losses_df, edge_nonvertical_df, edge_root_prok_df, edge_root_nonprok_df)

# Merge multiple mappings in ancestors to ensure consistent counts
edge_curr_to_update_agg <- c()
if (nrow(multiple_mappings) > 0) {
  edge_to_update <- edge_df %>% filter(origin_family %in% multiple_mappings$origin_family)
  edge_noupdate <- edge_df %>% filter(!origin_family %in% multiple_mappings$origin_family)
  
  for (i in 1:nrow(multiple_mappings)) {
    multiple_mapping <- multiple_mappings[i,]
    curr_node <- multiple_mapping$node
    curr_protein_ids <- unlist(strsplit(multiple_mapping$proteins, split=","))
    
    # Get ancestor nodes
    ancestor_nodes <- c()
    while (curr_node != (Ntip(tree)+1)) {
      curr_node <- get_ancestral_nodes(tree, curr_node, Nsplits=1) + Ntip(tree)
      ancestor_nodes <- c(ancestor_nodes, curr_node)
    }
    
    # Combine split parent PhROGs in ancestor nodes
    edge_curr_to_update <- edge_to_update %>% filter(origin_family %in% multiple_mapping$origin_family) %>% filter(child %in% ancestor_nodes & any(unlist(strsplit(protein_id, split=",")) %in% curr_protein_ids))
    edge_curr_no_update <- edge_to_update %>% filter(origin_family %in% multiple_mapping$origin_family) %>% filter(!(child %in% ancestor_nodes & any(unlist(strsplit(protein_id, split=",")) %in% curr_protein_ids)))
    edge_curr_to_update <- edge_curr_to_update %>% group_by(origin_family, child) %>% summarize(parent = max(parent), count = 1, mito_localization_prob = mean(mito_localization_prob), protein_id = paste0(sort(unlist(strsplit(protein_id, split=","))), collapse=","), is_singleton = 0, is_nonvertical = 0, .groups = "keep") %>% filter(!duplicated(protein_id))
    
    edge_to_update <- edge_to_update %>% filter(!origin_family %in% multiple_mapping$origin_family)
    edge_curr_to_update_agg <- rbind(edge_curr_to_update_agg, edge_curr_to_update)
    edge_to_update <- rbind(edge_to_update, edge_curr_to_update, edge_curr_no_update)
  }
  edge_df <- rbind(edge_noupdate, edge_to_update)
}

### Preprocess whole proteome data for gains and losses analysis
# Get absense data
homology_power_agg <- read.table(here("data/abSENSE_HMM", "absense_results.tsv"), sep="\t", header=TRUE)
colnames(homology_power_agg) <- c("PhROG_id", "human_accessions", "L", "R", "OG_mrca_label", "OG_outgroup_mrca_label", "outgroup_kingdom_species", "fraction_of_detectable_outgroup_kingdoms", "probability_of_detection_in_any_outgroup_species", "hmm_length", "r_squared", "comments")
homology_power_per_species <- read.table(here("data/abSENSE_HMM", "absense_results_per_species.tsv"), sep="\t", header=FALSE)
colnames(homology_power_per_species) <- c("PhROG_id", "taxid", "probability_of_detection")
homology_power_agg <- homology_power_agg %>% rowwise() %>% mutate(OG_id = gsub("_.*", "", PhROG_id))
# Mark OGs with Eukaryota as MRCA as 1
homology_power_agg$fraction_of_detectable_outgroup_kingdoms[homology_power_agg$OG_mrca_label == "Node34_Eukaryota"] <- 1
homology_power_agg$absense_label <- "Undetectable"
homology_power_agg$absense_label[homology_power_agg$fraction_of_detectable_outgroup_kingdoms > 0] <- "Detectable"
parent_progs_long$absense_label <- homology_power_agg$absense_label[match(parent_progs_long$PROG_id, homology_power_agg$PhROG_id)]
parent_progs_long_detectable_protein_ids <- unlist(strsplit(parent_progs_long$protein_id[parent_progs_long$absense_label == "Detectable"], ","))

## Generate gene table for all homologs (mito + nonmito). Include non-primary OGs
gene_table_all <- c()
for (i in 1:nrow(tax_levels_df)) {
  tax_levels_df_curr <- tax_levels_df[i,]
  print(tax_levels_df_curr$label)
  
  progs_long_filter <- progs_long_agg %>% filter(grepl(tax_levels_df_curr$label, PROG_id, fixed=TRUE)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))
  
  # Select PhROGs present at the tax level or inherited from earlier tax level
  if (tax_levels_df_curr$node_index != Ntip(tree)+1) {
    # If not the root, find rows that are present at the current tax level or were present in ancestor
    ancestor_node_index <- tree$edge[tree$edge[,2] == tax_levels_df_curr$node_index,1]
    gene_table_ancestor <- gene_table_all %>% filter(node == ancestor_node_index)
    
    parent_OG_protein_ids <- paste0(gene_table_ancestor$origin_family, "_", gene_table_ancestor$protein_id)
    OG_protein_ids <- paste0(progs_long_filter$OG_id, "_", progs_long_filter$protein_id)
    present_PhROG_ids <- unique(progs_long_filter$PROG_id[which(progs_long_filter$label == tax_levels_df_curr$label | OG_protein_ids %in% parent_OG_protein_ids)])
    progs_long_filter <- progs_long_filter %>% filter(PROG_id %in% present_PhROG_ids)
    
  } else {
    # If root, find rows that are present at the current tax level or Eukaryota_parent (to include duplications within LECA that have a post-LECA LCA)
    parent_OG_protein_ids <- paste0(parent_progs_long_eukaryota$OG_id, "_", parent_progs_long_eukaryota$protein_id)
    OG_protein_ids <- paste0(progs_long_filter$OG_id, "_", progs_long_filter$protein_id)
    present_PhROG_ids <- unique(progs_long_filter$PROG_id[which(progs_long_filter$label == tax_levels_df_curr$label | OG_protein_ids %in% parent_OG_protein_ids)])
    progs_long_filter <- progs_long_filter %>% filter(PROG_id %in% present_PhROG_ids)
    
    # Require that it be present in a euk supergroup.
    permitted_leca_duplication_LCA_nodelabels <- unique(c("Eukaryota", "Opimoda", "Amorphea_CRuMs", "Amorphea", "Obazoa", "Diphoda", "Diaphorectickes", "CAM_Haptista", "CAM"))
    progs_long_filter <- progs_long_filter %>% group_by(PROG_id) %>% filter(gsub("Node[0-9]+_", "", label) %in% permitted_leca_duplication_LCA_nodelabels)
    
    # Minimum species filter
    progs_long_filter <- progs_long_filter %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)
    
  }
  
  progs_long_filter <- progs_long_filter %>% rowwise() %>% mutate(mito_localization_prob_max = max(mito_localization_prob_mk, mito_localization_prob_parsimony, 0, na.rm=TRUE))
  
  gene_table_all_curr <- data.frame(node = tax_levels_df_curr$node_index, gene_id = progs_long_filter$PROG_id, origin_family = progs_long_filter$OG_id, protein_id = progs_long_filter$protein_id, mito_localization_prob = progs_long_filter$mito_localization_prob_max, stringsAsFactors = FALSE)
  gene_table_all <- rbind(gene_table_all, gene_table_all_curr)
}

# Get reference species whole proteomes. Include non-primary OGs
for (curr_taxid in eukaryote_reference_species_list) {
  ogs_long_curr <- ogs_long %>% filter(taxid == curr_taxid)
  deeploc_results_curr_species <- deeploc_results %>% filter(taxid == curr_taxid)
  deeploc_results_curr_species_mito_prob <- deeploc_results_curr_species$Mitochondrion[match(ogs_long_curr$accession, deeploc_results_curr_species$Protein_ID)]
  gene_table_curr <- data.frame(node = match(curr_taxid, tree_labels), gene_id = ogs_long_curr$accession, origin_family = ogs_long_curr$Orthogroup, protein_id = ogs_long_curr$accession, mito_localization_prob = deeploc_results_curr_species_mito_prob)
  gene_table_all <- rbind(gene_table_all, gene_table_curr)
}


# Add proteins from mtDNA-rich euks
gene_table_all <- rbind(gene_table_all, gene_table_add_mtdna_rich_euks_agg)
gene_table_all <- gene_table_all[!duplicated(gene_table_all),]


### Categorize losses
losses_agg <- c()
for (curr_label in tree_labels) {
  print(curr_label)
  selected_child_node <- which(tree_labels == curr_label)
  present_in_descendant <- edge_df %>% filter(child == selected_child_node)
  ancestor_index <- tree$edge[tree$edge[,2] == selected_child_node,1]
  if (length(ancestor_index) == 0) {
    next
  }
  present_in_ancestor <- edge_df %>% filter(child == ancestor_index)
  
  descendant_mito_protein_ids <- unlist(strsplit(present_in_descendant$protein_id, split=","))
  lost_from_mito <- present_in_ancestor %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% descendant_mito_protein_ids))
  
  if (selected_child_node > Ntip(tree)) {
    descendant_labels <- get_subtree_at_node(tree, selected_child_node - Ntip(tree))$subtree$tip.label
  } else {
    descendant_labels <- tree$tip.label[selected_child_node]
  }
  
  # Get proteins present in the descendant
  gene_table_all_descendant <- gene_table_all %>% filter(node == selected_child_node)
  gene_table_all_descendant_protein_ids <- unique(unlist(strsplit(gene_table_all_descendant$protein_id, ",")))
  
  # Complete loss: have no homologs in descendant
  complete_loss <- lost_from_mito %>% filter(!origin_family %in% gene_table_all_descendant$origin_family)
  complete_loss$type <- "Complete loss"
  complete_loss_protein_ids <- unlist(strsplit(complete_loss$protein_id, split=","))
  
  # Annotate complete losses using absense
  if (curr_label %in% species_tree$tip.label) {
    curr_clade_labels <- curr_label
  } else {
    subtree <- get_subtree_at_node(species_tree, curr_label)$subtree
    curr_clade_labels <- subtree$tip.label
  }
  parent_progs_long_curr <- parent_progs_long %>% filter(protein_id %in% complete_loss_protein_ids)
  homology_power_per_species_detectable <- homology_power_per_species %>% filter(PhROG_id %in% parent_progs_long_curr$PROG_id) %>% filter(taxid %in% curr_clade_labels) %>% filter(probability_of_detection > 0.95)
  PhROG_ids_any_detectable_homolog <- unique(homology_power_per_species_detectable$PhROG_id)
  parent_progs_long_curr_detectable_protein_ids <- parent_progs_long_curr$protein_id[parent_progs_long_curr$PROG_id %in% PhROG_ids_any_detectable_homolog]
  complete_loss_detectable <- complete_loss %>% rowwise() %>% filter(any(unlist(strsplit(protein_id, split=",")) %in% parent_progs_long_curr_detectable_protein_ids))
  complete_loss_undetectable <- complete_loss %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% parent_progs_long_curr_detectable_protein_ids))
  complete_loss_detectable$type <- "Complete loss of all homologs"
  complete_loss_undetectable$type <- "Possible loss (limited by homology detection)"

  # Lost all mito homologs: have only nonmito homologs in descendant and no direct descendants
  lost_all_mito_homologs <- lost_from_mito %>% rowwise() %>% filter(!origin_family %in% present_in_descendant$origin_family) %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% c(complete_loss_protein_ids, gene_table_all_descendant_protein_ids))) # MODIFIED
  lost_all_mito_homologs_protein_ids <- unlist(strsplit(lost_all_mito_homologs$protein_id, split=","))
  # Absense-HMM
  parent_progs_long_curr <- parent_progs_long %>% filter(protein_id %in% lost_all_mito_homologs_protein_ids)
  homology_power_per_species_detectable <- homology_power_per_species %>% filter(PhROG_id %in% parent_progs_long_curr$PROG_id) %>% filter(taxid %in% curr_clade_labels) %>% filter(probability_of_detection > 0.95)
  # Mark as detected if the parent PhROG has a detected homolog in the clade
  parent_PhROG_ids_curr <- unique(parent_progs_long_curr$PROG_id)
  parent_progs_long_has_homolog <- parent_progs_long %>% filter(PROG_id %in% parent_PhROG_ids_curr) %>% filter(taxid %in% curr_clade_labels)
  PhROG_ids_any_detectable_homolog <- unique(c(homology_power_per_species_detectable$PhROG_id, parent_progs_long_has_homolog$PROG_id))
  parent_progs_long_curr_detectable_protein_ids <- parent_progs_long_curr$protein_id[parent_progs_long_curr$PROG_id %in% PhROG_ids_any_detectable_homolog]
  lost_all_mito_homologs_detectable <- lost_all_mito_homologs %>% rowwise() %>% filter(any(unlist(strsplit(protein_id, split=",")) %in% parent_progs_long_curr_detectable_protein_ids))
  lost_all_mito_homologs_undetectable <- lost_all_mito_homologs %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% parent_progs_long_curr_detectable_protein_ids))
  lost_all_mito_homologs_detectable$type <- "Complete loss of mito homologs"
  lost_all_mito_homologs_undetectable$type <- "Possible loss (limited by homology detection)"

  # Relocated: have only nonmito homologs in descendant and has a direct descendant
  relocated <- lost_from_mito %>% rowwise() %>% filter(!origin_family %in% present_in_descendant$origin_family) %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% c(complete_loss_protein_ids, lost_all_mito_homologs_protein_ids))) # MODIFIED
  relocated_protein_ids <- unlist(strsplit(relocated$protein_id, split=","))
  relocated$type <- "Relocated to nonmito"

  # Remainder are partially lost from mito: retain a mito homolog
  partial_loss <- lost_from_mito %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% c(complete_loss_protein_ids, lost_all_mito_homologs_protein_ids, relocated_protein_ids)))
  partial_loss_protein_ids <- unlist(strsplit(partial_loss$protein_id, split=","))
  # Absense-HMM
  parent_progs_long_curr <- parent_progs_long %>% filter(protein_id %in% partial_loss_protein_ids)
  homology_power_per_species_detectable <- homology_power_per_species %>% filter(PhROG_id %in% parent_progs_long_curr$PROG_id) %>% filter(taxid %in% curr_clade_labels) %>% filter(probability_of_detection > 0.95)
  # Mark as detected if the parent PhROG has a detected homolog in the clade
  parent_PhROG_ids_curr <- unique(parent_progs_long_curr$PROG_id)
  parent_progs_long_has_homolog <- parent_progs_long %>% filter(PROG_id %in% parent_PhROG_ids_curr) %>% filter(taxid %in% curr_clade_labels)
  PhROG_ids_any_detectable_homolog <- unique(c(homology_power_per_species_detectable$PhROG_id, parent_progs_long_has_homolog$PROG_id))
  parent_progs_long_curr_detectable_protein_ids <- parent_progs_long_curr$protein_id[parent_progs_long_curr$PROG_id %in% PhROG_ids_any_detectable_homolog]
  partial_loss_detectable <- partial_loss %>% rowwise() %>% filter(any(unlist(strsplit(protein_id, split=",")) %in% parent_progs_long_curr_detectable_protein_ids))
  partial_loss_undetectable <- partial_loss %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% parent_progs_long_curr_detectable_protein_ids))
  partial_loss_detectable$type <- "Partial loss of mito homologs"
  partial_loss_undetectable$type <- "Possible loss (limited by homology detection)"
  
  losses_curr <- rbind(complete_loss_detectable, complete_loss_undetectable, lost_all_mito_homologs_detectable, lost_all_mito_homologs_undetectable, partial_loss_detectable, partial_loss_undetectable, relocated)
  losses_curr$label <- curr_label
  losses_agg <- rbind(losses_agg, losses_curr)
}

losses_summary <- losses_agg %>% group_by(label, type) %>% summarize(count = n())

losses_agg$type <- factor(losses_agg$type, levels = c("Possible loss (limited by homology detection)", "Complete loss of all homologs", "Complete loss of mito homologs", "Partial loss of mito homologs", "Relocated to nonmito"))
losses_summary$type <- factor(losses_summary$type, levels = c("Possible loss (limited by homology detection)", "Complete loss of all homologs", "Complete loss of mito homologs", "Partial loss of mito homologs", "Relocated to nonmito"))

losses_summary_species <- losses_summary %>% group_by(label) %>% summarize(count = sum(count))


### Categorize gains
hgt_gains_agg <- c()
retargeted_gains_agg <- c()
gene_birth_gains_detectable_agg <- c()
gene_birth_gains_notdetectable_agg <- c()
descendant_duplications_agg <- c()
duplication_counts <- c()

for (curr_label in tree_labels) {
  print(curr_label)
  selected_child_node <- which(tree_labels == curr_label)
  ancestor_index <- tree$edge[tree$edge[,2] == selected_child_node,1]
  if (length(ancestor_index) == 0) {
    next
  }
  present_in_ancestor <- edge_df %>% filter(child == ancestor_index)
  present_in_descendant <- edge_df %>% filter(child == selected_child_node)
  
  all_present_in_ancestor <- gene_table_all %>% filter(node == ancestor_index)
  all_present_in_ancestor_protein_ids <- unique(unlist(strsplit(all_present_in_ancestor$protein_id, ",")))
  
  # Find duplications
  descendant_duplications <- c()
  n_duplications <- 0
  for (i in 1:nrow(present_in_ancestor)) {
    present_in_ancestor_curr <- present_in_ancestor[i,]
    present_in_ancestor_curr_origin_family <- unique(present_in_ancestor_curr$origin_family)
    present_in_ancestor_curr_protein_ids <- unlist(strsplit(present_in_ancestor_curr$protein_id, ","))
    
    descendant_candidate_duplications <- present_in_descendant %>% rowwise() %>% filter(origin_family %in% present_in_ancestor_curr_origin_family) %>% filter(any(unlist(strsplit(protein_id, ",")) %in% present_in_ancestor_curr_protein_ids)) ## MODIFIED
    
    if (nrow(descendant_candidate_duplications) > 1) {
      descendant_duplications <- rbind(descendant_duplications, descendant_candidate_duplications)
      n_new_duplications <- nrow(descendant_candidate_duplications) - 1
      n_duplications <- n_duplications + n_new_duplications
    }
  }
  
  # Find all newly emerged mito proteins that are absent in ancestor
  gains_df <- edge_df %>% filter(parent == (selected_child_node + 0.25))
  # Exclude duplications from the other gains, if any overlap
  if (!is.null(descendant_duplications)) {
    descendant_duplications_protein_ids <- unlist(strsplit(descendant_duplications$protein_id, ","))
    gains_df <- gains_df %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, split=",")) %in% descendant_duplications_protein_ids))
  }
  
  gains_df_protein_ids <- unlist(strsplit(gains_df$protein_id, ","))
  
  # Find any existing ancestral proteins that were relocated to mito, using all phrogs (mito + nonmito) from above
  retargeted_gains <- gains_df %>% rowwise() %>% filter(any(unlist(strsplit(protein_id, ",")) %in% all_present_in_ancestor_protein_ids))
  retargeted_gains_protein_ids <- unlist(strsplit(retargeted_gains$protein_id, ","))
  
  # If there is an exact match to a prokaryote-origin Eukaryota_parent phrog, then indicates HGT
  parent_progs_long_prokorigin_curr <- parent_progs_long_prokorigin %>% filter(OG_id %in% gains_df$origin_family)
  parent_progs_long_prokorigin_hgt_gains <- parent_progs_long_prokorigin_curr %>% group_by(PROG_id) %>% filter(all(unlist(strsplit(protein_id, ",")) %in% gains_df_protein_ids))
  parent_progs_long_prokorigin_hgt_gains_protein_ids <- unlist(strsplit(parent_progs_long_prokorigin_hgt_gains$protein_id, ","))
  hgt_gains <- gains_df %>% rowwise() %>% filter(any(unlist(strsplit(protein_id, ",")) %in% parent_progs_long_prokorigin_hgt_gains_protein_ids)) %>% filter(!any(unlist(strsplit(protein_id, ",")) %in% retargeted_gains_protein_ids))
  hgt_gains_protein_ids <- unlist(strsplit(hgt_gains$protein_id, ","))
  
  # Remainder of gains are novel gene births. Split up into detectable or undetectable/uncertain using absense
  gene_birth_gains <- gains_df %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, ",")) %in% c(retargeted_gains_protein_ids, hgt_gains_protein_ids, descendant_duplications_protein_ids)))
  gene_birth_gains_detectable <- gene_birth_gains %>% rowwise() %>% filter(any(unlist(strsplit(protein_id, ",")) %in% parent_progs_long_detectable_protein_ids))
  gene_birth_gains_notdetectable <- gene_birth_gains %>% rowwise() %>% filter(!any(unlist(strsplit(protein_id, ",")) %in% parent_progs_long_detectable_protein_ids))
  
  if (!is.null(hgt_gains)) {
    hgt_gains$label <- curr_label
    hgt_gains_agg <- rbind(hgt_gains_agg, hgt_gains)
  }
  
  if (!is.null(retargeted_gains)) {
    retargeted_gains$label <- curr_label
    retargeted_gains_agg <- rbind(retargeted_gains_agg, retargeted_gains)
  }
  
  if (!is.null(gene_birth_gains_detectable)) {
    gene_birth_gains_detectable$label <- curr_label
    gene_birth_gains_detectable_agg <- rbind(gene_birth_gains_detectable_agg, gene_birth_gains_detectable)
  }
  
  if (!is.null(gene_birth_gains_notdetectable)) {
    gene_birth_gains_notdetectable$label <- curr_label
    gene_birth_gains_notdetectable_agg <- rbind(gene_birth_gains_notdetectable_agg, gene_birth_gains_notdetectable)
  }
  
  if (!is.null(descendant_duplications)) {
    descendant_duplications$label <- curr_label
    descendant_duplications_agg <- rbind(descendant_duplications_agg, descendant_duplications)
  }
  
  duplication_counts <- rbind(duplication_counts, data.frame(label = curr_label, count = n_duplications))
}

hgt_gains_agg$type <- "HGT"
retargeted_gains_agg$type <- "Relocated to mito"
gene_birth_gains_detectable_agg$type <- "Gene birth"
gene_birth_gains_notdetectable_agg$type <- "Possible gene birth (limited by homology detection)"
duplication_counts$type <- "Duplication of mito homolog"
descendant_duplications_agg$type <- "Duplication of mito homolog"

hgt_gains_summary <- hgt_gains_agg %>% group_by(label, type) %>% summarize(count = n())
retargeted_gains_summary <- retargeted_gains_agg %>% group_by(label, type) %>% summarize(count = n())
gene_birth_gains_detectable_summary <- gene_birth_gains_detectable_agg %>% group_by(label, type) %>% summarize(count = n())
gene_birth_gains_notdetectable_summary <- gene_birth_gains_notdetectable_agg %>% group_by(label, type) %>% summarize(count = n())

gains_summary <- rbind(hgt_gains_summary, retargeted_gains_summary, gene_birth_gains_detectable_summary, gene_birth_gains_notdetectable_summary, duplication_counts)

gains_summary$type <- factor(gains_summary$type, levels = c("Possible gene birth (limited by homology detection)", "Gene birth", "HGT", "Duplication of mito homolog", "Relocated to mito"))

gains_summary_species <- gains_summary %>% group_by(label) %>% summarize(label = unique(label), count = sum(count))


### Get ancestral gene counts
child_counts <- edge_df %>% filter(child %in% 1:length(tree_labels)) %>% group_by(child) %>% summarize(n_proteins = sum(count), n_OGs = length(unique(origin_family)))
child_counts$label <- tree_labels[child_counts$child]

# Get upper bound by absense
gains_summary_notdetectable <- gains_summary %>% filter(type == "Possible gene birth (limited by homology detection)")
losses_summary_notdetectable <- losses_summary %>% filter(type == "Possible loss (limited by homology detection)")

child_counts$possible_gains <- 0
child_counts$uncertain_losses <- 0

for (curr_label in tree$node.label) {
  if (curr_label == "Node34_Eukaryota") {
    # At the root, get the two daughter subtrees
    subtree <- get_subtree_at_node(tree, curr_label)$subtree
    subtree_labels <- c(subtree$tip.label, subtree$node.label)
    
    curr_node_index <- which(subtree$node.label == curr_label)
    daughter_nodes <- subtree$edge[subtree$edge[,1] == curr_node_index + Ntip(subtree),2]
    if (daughter_nodes[1] <= Ntip(subtree)) {
      subtree_1_labels <- subtree$tip.label[daughter_nodes[1]]
    } else {
      subtree_1 <- get_subtree_at_node(subtree, daughter_nodes[1] - Ntip(subtree))$subtree
      subtree_1_labels <- c(subtree_1$tip.label, subtree_1$node.label)
    }
    if (daughter_nodes[2] <= Ntip(subtree)) {
      subtree_2_labels <- subtree$tip.label[daughter_nodes[2]]
    } else {
      subtree_2 <- get_subtree_at_node(subtree, daughter_nodes[2] - Ntip(subtree))$subtree
      subtree_2_labels <- c(subtree_2$tip.label, subtree_2$node.label)
    }
  } else {
    # Bisect tree at current node, forming two subtrees
    subtree_1 <- get_subtree_at_node(tree, curr_label)$subtree
    subtree_1_labels <- c(subtree_1$tip.label, subtree_1$node.label)
    subtree_2 <- drop.tip(tree, subtree_1$tip.label)
    subtree_2_labels <- tree_labels[!tree_labels %in% subtree_1_labels]
    
    # Exclude current label
    subtree_1_labels <- subtree_1_labels[subtree_1_labels != curr_label]
    subtree_2_labels <- subtree_2_labels[subtree_2_labels != curr_label]
  }
  
  # Get undetectable gains
  gains_summary_notdetectable_descendants <- gains_summary_notdetectable %>% filter(label != curr_label)
  subtree1_gains_not_detectable <- sum(gains_summary_notdetectable_descendants$count[gains_summary_notdetectable_descendants$label %in% subtree_1_labels])
  subtree2_gains_not_detectable <- sum(gains_summary_notdetectable_descendants$count[gains_summary_notdetectable_descendants$label %in% subtree_2_labels])
  gains_not_detectable <- min(subtree1_gains_not_detectable, subtree2_gains_not_detectable) # a true gain should be an undetectable gain in both subtrees

  # Get undetectable losses
  losses_summary_notdetectable_curr <- losses_summary_notdetectable %>% filter(label == curr_label)
  losses_not_detectable <- sum(losses_summary_notdetectable_curr$count)
  
  child_counts$possible_gains[child_counts$label == curr_label] <- gains_not_detectable
  child_counts$uncertain_losses[child_counts$label == curr_label] <- losses_not_detectable
}

child_counts <- child_counts %>% rowwise() %>% mutate(possible_missing_proteins = possible_gains + uncertain_losses) %>% mutate(n_proteins_upper_bound = n_proteins + possible_missing_proteins)

## Get lower bound by parsimony on experimentally-defined mitoproteomes
mito_OG_ids <- unique(unlist(strsplit(gold_gene_accession_OG_id_df$OG_id[which(gold_gene_accession_OG_id_df$taxid %in% completed_mitoproteomes_species_list)], split=",")))

parsimony_ancestral_mito_copies_mat <- matrix(data = 0, nrow = length(mito_OG_ids), ncol = length(tree$node.label))
rownames(parsimony_ancestral_mito_copies_mat) <- mito_OG_ids
colnames(parsimony_ancestral_mito_copies_mat) <- tree$node.label
for (i in 1:length(mito_OG_ids)) {
  curr_OG_id <- mito_OG_ids[i]
  
  # MitoCarta only
  gold_gene_accession_OG_id_df_curr <- gold_gene_accession_OG_id_df %>% filter(OG_id == curr_OG_id) %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% group_by(taxid) %>% summarize(n_mito_copies = length(unique(gene_accession)))
  curr_copies_vec <- gold_gene_accession_OG_id_df_curr$n_mito_copies[match(tree$tip.label, gold_gene_accession_OG_id_df_curr$taxid)]
  present_species <- as.character(unique(gold_gene_accession_OG_id_df_curr$taxid))
  
  # Fill in zero for no homolog.
  curr_copies_vec[is.na(curr_copies_vec)] <- 0

  if (length(present_species) == 1) {
    # Gained at tip, so no ancestral copies. No update needed
    next
    
  } else if (length(unique(curr_copies_vec[!is.na(curr_copies_vec)])) == 1) {
    parsimony_ancestral_copies <- rep(mean(curr_copies_vec, na.rm=TRUE), ncol(parsimony_ancestral_mito_copies_mat))
    
  } else {
    curr_state_vec <- curr_copies_vec + 1 # add 1 to ensure states are positive integers
    
    ## Dollo parsimony with proportional transition matrix
    # Get subtree spanning the extant descendants
    mrca_index <- get_mrca_of_set(tree, present_species)
    gain_subtree <- get_subtree_at_node(tree, mrca_index - Ntip(tree))$subtree
    curr_state_vec <- curr_state_vec[match(gain_subtree$tip.label, tree$tip.label)]
    
    # Infer ancestral copies by parsimony model
    n_states <- max(curr_state_vec, na.rm=TRUE)
    cost_matrix <- matrix(data = NA, nrow=n_states, ncol=n_states)
    for (row in 1:n_states) {
      for (column in 1:n_states) {
        cost_matrix[row,column] <- abs(row - column)
      }
    }
    cost_matrix[1,2:n_states] <- Inf
    asr_result <- hsp_max_parsimony(gain_subtree, tip_states = curr_state_vec, transition_costs = cost_matrix)
    likelihoods_mat <- asr_result$likelihoods
    gain_subtree_node_indexes <- (Ntip(gain_subtree)+1):(Ntip(gain_subtree) + gain_subtree$Nnode)
    ancestral_likelihoods <- matrix(likelihoods_mat[gain_subtree_node_indexes,], nrow=length(gain_subtree_node_indexes)) # Get the state likelihoods at internal nodes
    parsimony_ancestral_copies_gain_subtree <- apply(ancestral_likelihoods, 1, function(row) sum(row * 0:(ncol(ancestral_likelihoods)-1)))
    parsimony_ancestral_copies <- rep(0, length(tree$node.label))
    parsimony_ancestral_copies[match(gain_subtree$node.label, tree$node.label)] <- parsimony_ancestral_copies_gain_subtree
  }
  
  parsimony_ancestral_mito_copies_mat[i,] <- parsimony_ancestral_copies
}

# Add proteins from mtDNA-rich euks 
gene_table_add_mtdna_rich_euks_agg_parsimony <- gene_table_add_mtdna_rich_euks_agg %>% mutate(label = gsub(".*_Node", "Node", gsub("_PhROG.*", "", gene_id))) %>% filter(!origin_family %in% mito_OG_ids)
parsimony_ancestral_copies_mtdna_added <- c()
for (curr_OG_id in unique(gene_table_add_mtdna_rich_euks_agg_parsimony$origin_family)) {
  gene_table_add_mtdna_rich_euks_agg_parsimony_curr <- gene_table_add_mtdna_rich_euks_agg_parsimony %>% filter(origin_family == curr_OG_id)
  parsimony_ancestral_copies <- as.numeric(colnames(parsimony_ancestral_mito_copies_mat) %in% gene_table_add_mtdna_rich_euks_agg_parsimony_curr$label)
  parsimony_ancestral_copies_mtdna_added <- rbind(parsimony_ancestral_copies_mtdna_added, parsimony_ancestral_copies)
}
parsimony_ancestral_mito_copies_mat <- rbind(parsimony_ancestral_mito_copies_mat, parsimony_ancestral_copies_mtdna_added)

parsimony_ancestral_mito_copies_mat_counts <- colSums(parsimony_ancestral_mito_copies_mat)
parsimony_ancestral_mito_copies_mat_counts_df <- data.frame(label = names(parsimony_ancestral_mito_copies_mat_counts), count = parsimony_ancestral_mito_copies_mat_counts)

child_counts$n_proteins_lower_bound <- parsimony_ancestral_mito_copies_mat_counts_df$count[match(child_counts$label, parsimony_ancestral_mito_copies_mat_counts_df$label)]

child_counts <- child_counts %>% mutate(n_proteins_CI = paste0(n_proteins, " (", round(n_proteins_lower_bound), "-", round(n_proteins_upper_bound), ")"))
child_counts$n_proteins_CI[is.na(child_counts$n_proteins_lower_bound)] <- child_counts$n_proteins_upper_bound[is.na(child_counts$n_proteins_lower_bound)]



