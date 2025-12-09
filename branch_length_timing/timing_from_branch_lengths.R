### Infer relative timings from branch lengths

# Load libraries
library(tidyverse)
library(ape)
library(castor)
library(ggtree)
library(ggplot2)
library(gplots)
library(dplyr)
library(RColorBrewer)
library(phytools)
theme_set(theme_classic())

### Read in datasets
# Read in experimental and mtDNA mito proteins
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2025.09.30.tsv"), sep="\t", header=TRUE)

## Read in human MitoPathways
# Read in mapping of human gene IDs to uniprot
human_id2training_mito <- read.table(here("data/orthogroups/idmapping", "human.id2training_uniprot_symbols.txt"), sep="\t", header=TRUE)
human_id2training_mito <- human_id2training_mito[order(human_id2training_mito$Symbol),]
# Read in human MitoPathways
mitopathways <- read.table(here("data/mito_orthogroups", "human.path2gene.txt"), sep="\t", header=TRUE)
mitopathways_human_id2training_mito <- merge(mitopathways, human_id2training_mito, by.x="Gene", by.y="Symbol")
mitopathways_human_id2training_mito$gene_accession <- paste0("9606_", mitopathways_human_id2training_mito$Entry)

mitopathways_OGs <- merge(gold_gene_accession_OG_id_df, mitopathways_human_id2training_mito, by="gene_accession")
mitopathways_OGs_summary <- mitopathways_OGs %>% group_by(OG_id, MitoPathwayHierarchy) %>% summarize(count = n(), .groups = "keep")
mitopathways_OGs_summary %>% filter(!grepl("\\.", MitoPathwayHierarchy))
mitopathways_OGs_summary_filter <- mitopathways_OGs_summary %>% filter(count == max(count))
mitopathways_OGs_summary_filter <- mitopathways_OGs_summary_filter[nrow(mitopathways_OGs_summary_filter):1,] # reverse order to get most specific labels
mitopathways_OGs_summary_filter <- mitopathways_OGs_summary_filter[!duplicated(mitopathways_OGs_summary_filter$OG_id),]
mitopathways_OGs_summary_all <- mitopathways_OGs %>% group_by(OG_id) %>% summarize(count = n(), MitoPathwayHierarchy_collapse = paste0(sort(unique(unlist(strsplit(MitoPathwayHierarchy, split="\\.")))), collapse="."), .groups = "keep")
mitopathways_OGs_summary_filter$MitoPathwayHierarchy_collapse <- mitopathways_OGs_summary_all$MitoPathwayHierarchy_collapse[match(mitopathways_OGs_summary_filter$OG_id, mitopathways_OGs_summary_all$OG_id)]

origin_table <- read.table(here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t")
colnames(origin_table) <- c("OG_id", "n_euk_species_largest_euk_clade", "n_prok_species_largest_prok_clade", "LCA_node", "LCA_node_euks", "origin_domain", "dropped_tips")
origin_table$OG_id <- gsub(".faa_clipkit.gappy.msa", "", origin_table$OG_id, fixed=TRUE)
prok_origin_OG_ids <- origin_table$OG_id[origin_table$origin_domain == "Prokaryote"]

# Read in OGs
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")

## Find pre-LECA duplications
completed_mitoproteomes_species_list <- c("9606", "559292", "3702", "1257118", "5689", "185431", "5741", "32595")
eukaryote_reference_species_list <- completed_mitoproteomes_species_list # default

## Set mito localization thresholds for Mk model and parsimony model
deeploc_thresholds <- read.csv(here("data/deeploc", "deeploc_thresholds.csv"))
deeploc_mito_threshold <- deeploc_thresholds$threshold[deeploc_thresholds$label == "Mitochondrion"]
mito_localization_prob_mk_threshold <- deeploc_mito_threshold
mito_localization_prob_parsimony_threshold <- 0.5

# Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

# Read in Eukaryota_parent PhROGs
phrogs_long_eukaryota_parent <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long_eukaryota_parent) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_parent <- phrogs_long_eukaryota_parent %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))

# Read in Eukaryota_child PhROGs
progs_long <- read.table(file.path("data/phylogenetically_resolved_orthogroups", "PhROGs_long", paste0("PhROGs_at_Node34_Eukaryota_long.tsv")), sep="\t", header=TRUE)
colnames(progs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota <- phrogs_long_eukaryota %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))

## Get reconstructed LECA whole proteome PhROGs/OGs
leca_OG_PhROG_ids <- read.table(here("data/reconstruction", "leca_PhROG_OG_ids_corespecies_deeplocmc.min5spp_Eukaryota.parent.2supergroups.min4spp.filter.txt"))$V1
phrogs_long_eukaryota_leca <- phrogs_long_eukaryota %>% filter(PROG_id %in% leca_OG_PhROG_ids)
phrogs_long_eukaryota_parent_leca <- phrogs_long_eukaryota_parent %>% filter(label == "Node34_Eukaryota")

# Require minimum number of species to support LECA nodes
n_minimum_species_leca <- 4
phrogs_long_eukaryota_leca <- phrogs_long_eukaryota_leca %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)
phrogs_long_eukaryota_parent_leca <- phrogs_long_eukaryota_parent_leca %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)

## Get reconstructed LECA mito PhROGs
leca_mito_PhROG_ids <- read.table(here("data/reconstruction", "leca_mito_PhROG_ids_corespecies_deeplocmc.min5spp_Eukaryota.parent.2supergroups.min4spp.filter.txt"))$V1
phrogs_long_eukaryota_leca_mito <- phrogs_long_eukaryota %>% filter(PROG_id %in% leca_mito_PhROG_ids)
phrogs_long_eukaryota_parent_leca_mito <- phrogs_long_eukaryota_parent_leca %>% group_by(PROG_id) %>% filter(mito_localization_prob_mk >= mito_localization_prob_mk_threshold | mito_localization_prob_parsimony >= mito_localization_prob_parsimony_threshold)
# Filter out DeepLocMC predictions below a minimum number of species with predicted mito protein
deeploc_results_exclude_OG_ids <- read.table(here('data/reconstruction', 'deeploc_results_exclude_lt5mitospp_OG_ids.txt'))$V1
phrogs_long_eukaryota_parent_leca_mito <- phrogs_long_eukaryota_parent_leca_mito %>% filter(!OG_id %in% deeploc_results_exclude_OG_ids)

### Identify duplications
# Get candidate duplications
phrogs_long_eukaryota_leca_dup <- phrogs_long_eukaryota_leca %>% group_by(OG_id) %>% summarize(n = length(unique(PROG_id))) %>% filter(n > 1)
phrogs_long_eukaryota_parent_leca_dup_candidate <- phrogs_long_eukaryota_parent_leca %>% filter(OG_id %in% phrogs_long_eukaryota_leca_dup$OG_id)

# Convert to wide form
phrogs_long_eukaryota_leca_wide <- phrogs_long_eukaryota_leca %>% group_by(PROG_id) %>% summarize(OG_id = unique(OG_id), protein_ids = paste0(protein_id, collapse=","), mito_localization_prob_mk = unique(mito_localization_prob_mk))
phrogs_long_eukaryota_parent_leca_dup_candidate_wide <- phrogs_long_eukaryota_parent_leca_dup_candidate %>% group_by(PROG_id) %>% summarize(OG_id = unique(OG_id), protein_ids = paste0(protein_id, collapse=","), mito_localization_prob_mk = unique(mito_localization_prob_mk))

find_parent <- function(query_protein_ids, target_protein_ids) {
  query_protein_ids <- unique(unlist(strsplit(query_protein_ids, split=",")))
  shared_protein_ids <- intersect(target_protein_ids, query_protein_ids)
  return(length(shared_protein_ids) > 0)
}

# Iterate through parent PhROGs to find duplicated child PhROGs
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$child_protein_id <- ""
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$child_mito_localization_prob_mk <- NA
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$n_child_phrogs <- NA
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$n_child_mito_phrogs <- NA
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$n_child_nonmito_phrogs <- NA
for (i in 1:nrow(phrogs_long_eukaryota_parent_leca_dup_candidate_wide)) {
  # Filter for matching OG
  phrogs_long_eukaryota_leca_wide_curr <- phrogs_long_eukaryota_leca_wide %>% filter(OG_id %in% phrogs_long_eukaryota_parent_leca_dup_candidate_wide$OG_id[i])
  
  # Find shared protein ids in parent and child
  target_protein_ids <- unique(unlist(strsplit(phrogs_long_eukaryota_parent_leca_dup_candidate_wide$protein_ids[i], split=",")))
  parent_index <- which(sapply(phrogs_long_eukaryota_leca_wide_curr$protein_ids, find_parent, target_protein_ids))
  
  # Update
  phrogs_long_eukaryota_parent_leca_dup_candidate_wide$child_protein_id[i] <- paste0(phrogs_long_eukaryota_leca_wide_curr$protein_ids[parent_index], collapse="+")
  phrogs_long_eukaryota_parent_leca_dup_candidate_wide$child_mito_localization_prob_mk[i] <- paste0(phrogs_long_eukaryota_leca_wide_curr$mito_localization_prob_mk[parent_index], collapse="+")
  phrogs_long_eukaryota_parent_leca_dup_candidate_wide$n_child_phrogs[i] <- length(parent_index)
  phrogs_long_eukaryota_parent_leca_dup_candidate_wide$n_child_mito_phrogs[i] <- sum(phrogs_long_eukaryota_leca_wide_curr$PROG_id[parent_index] %in% phrogs_long_eukaryota_leca_mito$PROG_id)
  phrogs_long_eukaryota_parent_leca_dup_candidate_wide$n_child_nonmito_phrogs[i] <- sum(!phrogs_long_eukaryota_leca_wide_curr$PROG_id[parent_index] %in% phrogs_long_eukaryota_leca_mito$PROG_id)
}

# Assign mito localization states to candidate duplication parents
phrogs_long_eukaryota_parent_leca_dup_candidate_wide <- phrogs_long_eukaryota_parent_leca_dup_candidate_wide %>% mutate(parent_is_mito = PROG_id %in% phrogs_long_eukaryota_parent_leca_mito$PROG_id)

## Assign GO cellular components to parents. Assign GO term if >50% of annotated proteins have the GO term.
go_annot_orthogroups <- readRDS(here('data/annotation/GO', 'GO_cellular_components_refined_OG_mapping.rds'))
phrogs_long_eukaryota_parent_leca_go <- merge(phrogs_long_eukaryota_parent_leca, go_annot_orthogroups, by.x=c("protein_id", "OG_id"), by.y=c("accessions", "OG_id"))
# Map children to parent GO terms
go_name_to_id_table <- read.csv(here('data/annotation/GO', 'GO_cellular_components_name_to_id_children.csv'))
phrogs_long_eukaryota_parent_leca_go$goId_parent <- ""
for (i in 1:nrow(go_name_to_id_table)) {
  goIds_curr <- c(go_name_to_id_table$goId_parent[i], unlist(strsplit(go_name_to_id_table$goId_children[i], split=",")))
  phrogs_long_eukaryota_parent_leca_go$goId_parent[which(phrogs_long_eukaryota_parent_leca_go$goId %in% goIds_curr)] <- paste0(phrogs_long_eukaryota_parent_leca_go$goId_parent[which(phrogs_long_eukaryota_parent_leca_go$goId %in% goIds_curr)], ",", go_name_to_id_table$goId_parent[i])
}
phrogs_long_eukaryota_parent_leca_go <- phrogs_long_eukaryota_parent_leca_go %>% separate_rows(goId_parent, sep=",")
phrogs_long_eukaryota_parent_leca_go <- phrogs_long_eukaryota_parent_leca_go %>% filter(goId_parent != "")
# Remove duplicated protein-GO pairs
phrogs_long_eukaryota_parent_leca_go <- phrogs_long_eukaryota_parent_leca_go[!duplicated(cbind(phrogs_long_eukaryota_parent_leca_go$protein_id, phrogs_long_eukaryota_parent_leca_go$goId_parent)),]
phrogs_long_eukaryota_parent_leca_go_counts <- phrogs_long_eukaryota_parent_leca_go %>% group_by(PROG_id) %>% mutate(n_proteins_total = length(unique(protein_id)), n_GO_proteins = length(unique(protein_id)), n_species = length(unique(taxonId))) %>% group_by(PROG_id, taxonId) %>% mutate(weight = 1 / length(unique(protein_id))) %>% group_by(PROG_id, goId_parent) %>% summarize(n_proteins_total = unique(n_proteins_total), n_GO_proteins = unique(n_GO_proteins), n_species = unique(n_species), sum_weight = sum(weight)) %>% mutate(normalized_weight = sum_weight / n_species)
# Set threshold to assign GO term
phrogs_long_eukaryota_parent_leca_go_counts <- phrogs_long_eukaryota_parent_leca_go_counts %>% filter(normalized_weight > 0.5)

phrogs_long_eukaryota_parent_leca_go_summary <- phrogs_long_eukaryota_parent_leca_go_counts %>% group_by(PROG_id) %>% summarize(goIds = paste0(goId_parent, collapse=","))
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$goIds <- NA
phrogs_long_eukaryota_parent_leca_dup_candidate_wide$goIds[phrogs_long_eukaryota_parent_leca_dup_candidate_wide$PROG_id %in% phrogs_long_eukaryota_parent_leca_go_summary$PROG_id] <- phrogs_long_eukaryota_parent_leca_go_summary$goIds[match(phrogs_long_eukaryota_parent_leca_dup_candidate_wide$PROG_id[phrogs_long_eukaryota_parent_leca_dup_candidate_wide$PROG_id %in% phrogs_long_eukaryota_parent_leca_go_summary$PROG_id], phrogs_long_eukaryota_parent_leca_go_summary$PROG_id)]
# Add columns for selected GO terms (parent id + child ids)
for (i in 1:nrow(go_name_to_id_table)) {
  goIds_curr <- go_name_to_id_table$goId_parent[i]
  phrogs_long_eukaryota_parent_leca_dup_candidate_wide[go_name_to_id_table$name[i]] <- sapply(phrogs_long_eukaryota_parent_leca_dup_candidate_wide$goIds, function(s, goIds_curr) {
    any(unlist(strsplit(s, split=",")) %in% goIds_curr)
  }, goIds_curr)
}

## Assign GO cellular components to children. Assign GO term if >=50% of annotated proteins have the GO term.
phrogs_long_eukaryota_leca_go <- merge(phrogs_long_eukaryota_leca, go_annot_orthogroups, by.x=c("protein_id", "OG_id"), by.y=c("accessions", "OG_id"))
# Map children to parent GO terms
phrogs_long_eukaryota_leca_go$goId_parent <- ""
for (i in 1:nrow(go_name_to_id_table)) {
  goIds_curr <- c(go_name_to_id_table$goId_parent[i], unlist(strsplit(go_name_to_id_table$goId_children[i], split=",")))
  phrogs_long_eukaryota_leca_go$goId_parent[which(phrogs_long_eukaryota_leca_go$goId %in% goIds_curr)] <- paste0(phrogs_long_eukaryota_leca_go$goId_parent[which(phrogs_long_eukaryota_leca_go$goId %in% goIds_curr)], ",", go_name_to_id_table$goId_parent[i])
}
phrogs_long_eukaryota_leca_go <- phrogs_long_eukaryota_leca_go %>% separate_rows(goId_parent, sep=",")
phrogs_long_eukaryota_leca_go <- phrogs_long_eukaryota_leca_go %>% filter(goId_parent != "")

# Remove duplicated protein-GO pairs
phrogs_long_eukaryota_leca_go <- phrogs_long_eukaryota_leca_go[!duplicated(cbind(phrogs_long_eukaryota_leca_go$protein_id, phrogs_long_eukaryota_leca_go$goId_parent)),]
phrogs_long_eukaryota_leca_go_counts <- phrogs_long_eukaryota_leca_go %>% group_by(PROG_id) %>% mutate(n_proteins_total = length(unique(protein_id)), n_GO_proteins = length(unique(protein_id)), n_species = length(unique(taxonId))) %>% group_by(PROG_id, taxonId) %>% mutate(weight = 1 / length(unique(protein_id))) %>% group_by(PROG_id, goId_parent) %>% summarize(n_proteins_total = unique(n_proteins_total), n_GO_proteins = unique(n_GO_proteins), n_species = unique(n_species), sum_weight = sum(weight)) %>% mutate(normalized_weight = sum_weight / n_species)

# Filter by >= 0.5 normalized weight per species to assign GO term
phrogs_long_eukaryota_leca_go_counts <- phrogs_long_eukaryota_leca_go_counts %>% filter(normalized_weight >= 0.5)
phrogs_long_eukaryota_leca_go_summary <- phrogs_long_eukaryota_leca_go_counts %>% group_by(PROG_id) %>% summarize(goIds = paste0(goId_parent, collapse=","))

# Mark parents with multiple children as duplicated
phrogs_long_eukaryota_parent_leca_dup <- phrogs_long_eukaryota_parent_leca_dup_candidate_wide %>% filter(n_child_phrogs > 1)

## Get the set of LECA OGs with prokaryote origin or pre-LECA duplication that can be analyzed
phrogs_long_eukaryota_parent_leca_dup_OG_ids <- unique(phrogs_long_eukaryota_parent_leca_dup$OG_id)
phrogs_long_eukaryota_parent_leca_OG_ids <- unique(phrogs_long_eukaryota_parent_leca$OG_id)
phrogs_long_eukaryota_parent_leca_prok_OG_ids <- phrogs_long_eukaryota_parent_leca_OG_ids[phrogs_long_eukaryota_parent_leca_OG_ids %in% prok_origin_OG_ids]
phrogs_long_eukaryota_parent_leca_prok_or_dup_OG_ids <- unique(c(phrogs_long_eukaryota_parent_leca_prok_OG_ids, phrogs_long_eukaryota_parent_leca_dup_OG_ids))

### LECA branch length timings
# Branch length analysis
library(ggridges)

bl_df <- read.table(here("data/branch_length_timing", "branch_lengths_Node34_Eukaryota_combined.tsv"), sep="\t", header=FALSE)
colnames(bl_df) <- c("OG_id", "parent_PhROG_id", "child_PhROG_id", "median_euk_clade_bl", "mean_euk_clade_bl", "sd_euk_clade_bl", "stem_bl", "normalized_stem_bl", "duplication_bl", "normalized_duplication_bl", "oldest_duplication_bl", "normalized_oldest_duplication_bl")

# Filter for LECA nodes (strict! require that clade LCA is LECA to guarantee existence of a LECA node and has minimum number of species)
Node34_Eukaryota_PhROG_ids <- unique(phrogs_long_eukaryota$PROG_id[phrogs_long_eukaryota$label == "Node34_Eukaryota"])
bl_df <- bl_df %>% filter(child_PhROG_id %in% Node34_Eukaryota_PhROG_ids)

# Filter for LECA OGs that have prokaryote origin or pre-LECA duplication 
bl_df <- bl_df %>% filter(OG_id %in% phrogs_long_eukaryota_parent_leca_prok_or_dup_OG_ids)

# Label with mito localization
bl_df <- bl_df %>% mutate(is_MitoCarta = child_PhROG_id %in% leca_mito_PhROG_ids)

# Log normalize
bl_df <- bl_df %>% mutate(neglog10_normalized_stem_bl = -log10(normalized_stem_bl), neglog10_normalized_duplication_bl = -log10(normalized_duplication_bl), neglog10_normalized_oldest_duplication_bl = -log10(normalized_oldest_duplication_bl))

## Combined length by GO cellular compartment
bl_stem <- bl_df %>% filter(!is.na(neglog10_normalized_stem_bl) & !is.infinite(neglog10_normalized_stem_bl))
bl_dup <- bl_df %>% filter(!is.na(neglog10_normalized_oldest_duplication_bl) & !is.infinite(neglog10_normalized_oldest_duplication_bl))
bl_dup <- bl_dup %>% filter(is.na(neglog10_normalized_stem_bl))
bl_combined <- data.frame(OG_id = c(bl_stem$OG_id, bl_dup$OG_id), parent_PhROG_id = c(bl_stem$parent_PhROG_id, bl_dup$parent_PhROG_id), child_PhROG_id = c(bl_stem$child_PhROG_id, bl_dup$child_PhROG_id), type = c(rep("Acquisition", nrow(bl_stem)), rep("Duplication", nrow(bl_dup))), is_MitoCarta = c(bl_stem$is_MitoCarta, bl_dup$is_MitoCarta), neglog10_normalized_branch_length = c(bl_stem$neglog10_normalized_stem_bl, bl_dup$neglog10_normalized_oldest_duplication_bl))
bl_combined$type <- factor(bl_combined$type, levels = c("Acquisition", "Duplication"))
bl_go <- merge(bl_combined, phrogs_long_eukaryota_leca_go_counts, by.x="child_PhROG_id", by.y="PROG_id")
bl_go$GO_name <- go_name_to_id_table$name[match(bl_go$goId_parent, go_name_to_id_table$goId_parent)]
temp <- bl_go %>% filter(!is.na(neglog10_normalized_branch_length) & !is.infinite(neglog10_normalized_branch_length))
temp <- temp %>% filter(!GO_name %in% c("plastid", "cell wall", "mitochondrion", "mitochondrial ribosome"))
bl_MitoCarta <- bl_combined %>% filter(is_MitoCarta) %>% filter(!is.na(neglog10_normalized_branch_length) & !is.infinite(neglog10_normalized_branch_length))
bl_MitoCarta$GO_name <- "mitochondrion"
temp <- rbind(temp[,colnames(temp) %in% colnames(bl_MitoCarta)], bl_MitoCarta)
temp <- temp %>% group_by(GO_name) %>% filter(!duplicated(child_PhROG_id))
temp_median <- temp %>% group_by(GO_name) %>% summarize(median = median(neglog10_normalized_branch_length))
temp_median <- temp_median[order(temp_median$median, decreasing=TRUE),]
temp$GO_name <- factor(temp$GO_name, levels = c(temp_median$GO_name))
temp_count <- temp %>% group_by(GO_name) %>% summarize(n = n())

# All-vs-all pairwise Wilcoxon rank-sum test, BH correction
pairwise.wilcox.test(temp$neglog10_normalized_branch_length, temp$GO_name, p.adjust.method = "BH", exact = FALSE)

# pdf('combined_length_GO_cellular_compartments_LECA.pdf', height = 4, width = 4)
# ggplot(temp, aes(x = neglog10_normalized_branch_length, y = GO_name, fill = GO_name)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") + geom_text(
#     data = temp_count,
#     aes(x = Inf, label = n),
#     hjust = 1, # Align text to the right
#     nudge_y = 0.25, # Move the labels slightly above the ridgeline
#     size = 2,
#     color = "black" 
#   ) + xlim(-1,3)
# dev.off()

## Stem length by prokaryotic donor
combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long <- read.table(here("data/horizontal_gene_transfer", "HGT_Node34_Eukaryota_long.tsv"), sep="\t", header=TRUE)
colnames(combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long) <- c("OG_id", "PhROG_id", "count", "fraction_euk_species", "sister_clade_fraction_euk_species", "cousin_clade_fraction_euk_species", "HGT_sister", "HGT_self_sister", "taxon", "bool_in_mito")
temp <- bl_df %>% filter(!is.na(neglog10_normalized_stem_bl) & !is.infinite(neglog10_normalized_stem_bl))
temp$prokaryotic_donor <- as.character(combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long$taxon[match(temp$parent_PhROG_id, combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long$PhROG_id)])
temp$prokaryotic_donor[is.na(temp$prokaryotic_donor)] <- "Indeterminate"
temp <- temp %>% filter(prokaryotic_donor != "Indeterminate")
other_bacteria_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Bacteria"]
other_bacteria_ids <- other_bacteria_ids[other_bacteria_ids != "Alphaproteobacteria"]
other_archaea_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Archaea"]
other_archaea_ids <- other_archaea_ids[other_archaea_ids != "Asgardgroup"]
temp$prokaryotic_donor[temp$prokaryotic_donor %in% other_bacteria_ids] <- "Other Bacteria"
temp$prokaryotic_donor[temp$prokaryotic_donor %in% other_archaea_ids] <- "Other Archaea"
temp_median <- temp %>% group_by(prokaryotic_donor) %>% summarize(median = median(neglog10_normalized_stem_bl))
temp_median <- temp_median[order(temp_median$median, decreasing=TRUE),]
temp$prokaryotic_donor <- factor(temp$prokaryotic_donor, levels = c(temp_median$prokaryotic_donor))
temp_count <- temp %>% group_by(prokaryotic_donor) %>% summarize(n = n())
stem_bl_alpha <- temp %>% filter(prokaryotic_donor == "Alphaproteobacteria")

# All-vs-all pairwise Wilcoxon rank-sum test, BH correction
pairwise.wilcox.test(temp$neglog10_normalized_stem_bl, temp$prokaryotic_donor, p.adjust.method = "BH", exact = FALSE)

# pdf('stem_length_prokaryotic_donor_LECA.pdf', height = 6, width = 6)
# ggplot(temp, aes(x = neglog10_normalized_stem_bl, y = prokaryotic_donor, fill = prokaryotic_donor)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   geom_text(
#     data = temp_count,
#     aes(x = Inf, label = n),
#     hjust = 1, # Align text to the right
#     nudge_y = 0.25, # Move the labels slightly above the ridgeline
#     size = 2,
#     color = "black" 
#   ) + xlim(-1,3)
# dev.off()


## Combined length by MitoPathway per gene
bl_stem <- bl_df %>% filter(!is.na(neglog10_normalized_stem_bl) & !is.infinite(neglog10_normalized_stem_bl))
bl_dup <- bl_df %>% filter(!is.na(neglog10_normalized_oldest_duplication_bl) & !is.infinite(neglog10_normalized_oldest_duplication_bl))
bl_dup <- bl_dup %>% filter(is.na(neglog10_normalized_stem_bl))
bl_combined <- data.frame(OG_id = c(bl_stem$OG_id, bl_dup$OG_id), parent_PhROG_id = c(bl_stem$parent_PhROG_id, bl_dup$parent_PhROG_id), child_PhROG_id = c(bl_stem$child_PhROG_id, bl_dup$child_PhROG_id), type = c(rep("Acquisition", nrow(bl_stem)), rep("Duplication", nrow(bl_dup))), is_MitoCarta = c(bl_stem$is_MitoCarta, bl_dup$is_MitoCarta), neglog10_normalized_branch_length = c(bl_stem$neglog10_normalized_stem_bl, bl_dup$neglog10_normalized_oldest_duplication_bl))
bl_combined$type <- factor(bl_combined$type, levels = c("Acquisition", "Duplication"))
bl_combined <- bl_combined %>% filter(is_MitoCarta)

# Get human mitopathways
mitopathways_reduce <- mitopathways_human_id2training_mito[,c("Gene", "MitoPathway", "gene_accession")]
phrogs_long_eukaryota_leca_mitopathways <- merge(phrogs_long_eukaryota_leca, mitopathways_reduce, by.x="protein_id", by.y="gene_accession", all.x = TRUE)
bl_combined_mitopathway <- merge(bl_combined, phrogs_long_eukaryota_leca_mitopathways, by.x=c("OG_id", "child_PhROG_id"), by.y=c("OG_id", "PROG_id"))

# Impute MitoPathways for unlabeled mito nodes
fraction_children_labeled_mitopathway_threshold <- 1
bl_combined_mitopathway_parent_impute <- bl_combined_mitopathway %>% filter(!is.na(MitoPathway)) %>% group_by(parent_PhROG_id) %>% mutate(n_children_labeled = length(unique(child_PhROG_id))) %>% group_by(parent_PhROG_id, MitoPathway) %>% mutate(n_children_mitopathway = length(unique(child_PhROG_id))) %>% mutate(fraction_children_labeled_mitopathway = n_children_mitopathway / n_children_labeled) %>% filter(fraction_children_labeled_mitopathway >= fraction_children_labeled_mitopathway_threshold) %>% group_by(parent_PhROG_id) %>% summarize(MitoPathway = paste0(MitoPathway, collapse=","))
bl_combined_mitopathway_labeled <- bl_combined_mitopathway %>% filter(!is.na(MitoPathway))
bl_combined_mitopathway_missing <- bl_combined_mitopathway %>% filter(is.na(MitoPathway))
bl_combined_mitopathway_missing$MitoPathway <- bl_combined_mitopathway_parent_impute$MitoPathway[match(bl_combined_mitopathway_missing$parent_PhROG_id, bl_combined_mitopathway_parent_impute$parent_PhROG_id)]
bl_combined_mitopathway_missing <- bl_combined_mitopathway_missing %>% filter(!is.na(MitoPathway))
bl_combined_mitopathway_missing <- bl_combined_mitopathway_missing %>% separate_rows(MitoPathway, sep=",")
bl_combined_mitopathway_missing$Gene[gsub("_.*", "", bl_combined_mitopathway_missing$protein_id) %in% eukaryote_reference_species_list] <- bl_combined_mitopathway_missing$protein_id[gsub("_.*", "", bl_combined_mitopathway_missing$protein_id) %in% eukaryote_reference_species_list]
bl_combined_mitopathway_labeled$mitopathway_type <- "Human_MitoPathway"
bl_combined_mitopathway_missing$mitopathway_type <- "Imputed_MitoPathway"
bl_combined_mitopathway <- rbind(bl_combined_mitopathway_labeled, bl_combined_mitopathway_missing)

bl_combined_mitopathway <- bl_combined_mitopathway %>% filter(!is.na(MitoPathway)) %>% group_by(parent_PhROG_id, child_PhROG_id, MitoPathway, type, neglog10_normalized_branch_length) %>% summarize(Gene = paste0(sort(unique(Gene)), collapse=", "))

# Add all mito proteins
bl_combined_allmito <- bl_combined
bl_combined_allmito$MitoPathway <- "mitochondrion"
bl_combined_allmito$Gene <- bl_combined_allmito$child_PhROG_id
bl_combined_mitopathway <- rbind(bl_combined_mitopathway, bl_combined_allmito[,colnames(bl_combined_allmito) %in% colnames(bl_combined_mitopathway)])

bl_combined_mitopathway_select <- bl_combined_mitopathway
bl_combined_mitopathway_select <- bl_combined_mitopathway_select[order(bl_combined_mitopathway_select$neglog10_normalized_branch_length, decreasing=FALSE),]
selected_mitopathways <- c("Protein_import_sorting_and_homeostasis", "Small_molecule_transport", "OXPHOS", "Mitochondrial_central_dogma", "Metabolism")
bl_combined_mitopathway_select <- bl_combined_mitopathway_select %>% filter(MitoPathway %in% selected_mitopathways)
bl_combined_mitopathway_select <- bl_combined_mitopathway_select %>% group_by(MitoPathway) %>% filter(!duplicated(Gene))
temp_median <- bl_combined_mitopathway_select %>% group_by(MitoPathway) %>% summarize(median = median(neglog10_normalized_branch_length))
temp_median <- temp_median[order(temp_median$median, decreasing=TRUE),]
bl_combined_mitopathway_select$MitoPathway <- factor(bl_combined_mitopathway_select$MitoPathway, levels = c(temp_median$MitoPathway))
temp_count <- bl_combined_mitopathway_select %>% group_by(MitoPathway) %>% summarize(n = n())

# All-vs-all pairwise Wilcoxon rank-sum test, BH correction
pairwise.wilcox.test(bl_combined_mitopathway_select$neglog10_normalized_branch_length, bl_combined_mitopathway_select$MitoPathway, p.adjust.method = "BH", exact = FALSE)

# pdf('combined_length_MitoPathways_human_imputed.threshold1_LECA.pdf', height = 6, width = 8)
# ggplot(data = bl_combined_mitopathway_select, aes(x = neglog10_normalized_branch_length, y = MitoPathway, fill = MitoPathway)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   geom_text(
#     data = temp_count,
#     aes(x = Inf, label = n),
#     hjust = 1, # Align text to the right
#     nudge_y = 0.25, # Move the labels slightly above the ridgeline
#     size = 2,
#     color = "black" 
#   ) + xlim(-1,3)
# dev.off()


### Archaeplastida branch length timings

## Get OGs with Archaeplastida nodes for branch length analysis
phrogs_long_archaeplastida <- read.table(file.path("data/phylogenetically_resolved_orthogroups", "PhROGs_long", paste0("PhROGs_at_Node39_Archaeplastida_long.tsv")), sep="\t", header=TRUE)
colnames(phrogs_long_archaeplastida) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_archaeplastida <- phrogs_long_archaeplastida %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))
phrogs_long_archaeplastida_ancestor <- phrogs_long_archaeplastida %>% filter(label == "Node39_Archaeplastida")

bl_df <- read.table(here("data/branch_length_timing", "branch_lengths_Node39_Archaeplastida_combined.tsv"), sep="\t", header=FALSE)
colnames(bl_df) <- c("OG_id", "parent_PhROG_id", "child_PhROG_id", "median_euk_clade_bl", "mean_euk_clade_bl", "sd_euk_clade_bl", "stem_bl", "normalized_stem_bl", "duplication_bl", "normalized_duplication_bl", "oldest_duplication_bl", "normalized_oldest_duplication_bl")

# Filter for Archaeplastida nodes (strict! require that clade LCA is Archaeplastida to guarantee existence of an Archaeplastida node)
bl_df <- bl_df %>% filter(child_PhROG_id %in% phrogs_long_archaeplastida_ancestor$PROG_id)

## Get mito PhROGs
# Filter by mito localization prob by Mk or parsimony
phrogs_long_archaeplastida_mito <- phrogs_long_archaeplastida %>% group_by(PROG_id) %>% filter(mito_localization_prob_mk >= mito_localization_prob_mk_threshold | mito_localization_prob_parsimony >= mito_localization_prob_parsimony_threshold)
bl_df <- bl_df %>% mutate(is_MitoCarta = child_PhROG_id %in% phrogs_long_archaeplastida_mito$PROG_id)

bl_df <- bl_df %>% mutate(neglog10_normalized_stem_bl = -log10(normalized_stem_bl), neglog10_normalized_duplication_bl = -log10(normalized_duplication_bl), neglog10_normalized_oldest_duplication_bl = -log10(normalized_oldest_duplication_bl))

## Stem length by prokaryotic donor
temp <- bl_df %>% filter(!is.na(neglog10_normalized_stem_bl) & !is.infinite(neglog10_normalized_stem_bl))
combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long <- read.table(here("data/horizontal_gene_transfer", "HGT_Node39_Archaeplastida_long.tsv"), sep="\t", header=TRUE)
colnames(combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long) <- c("OG_id", "PhROG_id", "count", "fraction_euk_species", "sister_clade_fraction_euk_species", "cousin_clade_fraction_euk_species", "HGT_sister", "HGT_self_sister", "taxon", "bool_in_mito")
temp$prokaryotic_donor <- as.character(combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long$taxon[match(temp$parent_PhROG_id, combined_hgt_phrog_HGT_self_sister_leca_prok_origin_all_long$PhROG_id)])
temp$prokaryotic_donor[is.na(temp$prokaryotic_donor)] <- "Indeterminate"
temp <- temp %>% filter(prokaryotic_donor != "Indeterminate")
other_bacteria_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Bacteria"]
other_bacteria_ids <- other_bacteria_ids[!other_bacteria_ids %in% c("Alphaproteobacteria", "CyanobacteriotaMelainabacteriagroup")]
other_archaea_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Archaea"]
other_archaea_ids <- other_archaea_ids[other_archaea_ids != "Asgardgroup"]
temp$prokaryotic_donor[temp$prokaryotic_donor %in% other_bacteria_ids] <- "Other Bacteria"
temp$prokaryotic_donor[temp$prokaryotic_donor %in% other_archaea_ids] <- "Other Archaea"
temp_median <- temp %>% group_by(prokaryotic_donor) %>% summarize(median = median(neglog10_normalized_stem_bl))
temp_median <- temp_median[order(temp_median$median, decreasing=TRUE),]
temp$prokaryotic_donor <- factor(temp$prokaryotic_donor, levels = c(temp_median$prokaryotic_donor))
temp_count <- temp %>% group_by(prokaryotic_donor) %>% summarize(n = n())

# All-vs-all pairwise Wilcoxon rank-sum test, BH correction
pairwise.wilcox.test(temp$normalized_stem_bl, temp$prokaryotic_donor, p.adjust.method = "BH", exact = FALSE)

# pdf('stem_length_prokaryotic_donor_Archaeplastida.pdf', height = 6, width = 6)
# ggplot(temp, aes(x = neglog10_normalized_stem_bl, y = prokaryotic_donor, fill = prokaryotic_donor)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   geom_text(
#     data = temp_count,
#     aes(x = Inf, label = n),
#     hjust = 1, # Align text to the right
#     nudge_y = 0.25, # Move the labels slightly above the ridgeline
#     size = 2,
#     color = "black" 
#   ) + xlim(-1.5,3)
# dev.off()


## Combined length by GO cellular compartment
# Assign GO cellular components to children. Assign GO term if >50% of annotated proteins have the GO term.
go_annot_orthogroups <- readRDS(here("data/annotation/GO", "GO_cellular_components_Archaeplastida_refined_OG_mapping.rds"))
phrogs_long_archaeplastida_go <- merge(phrogs_long_archaeplastida, go_annot_orthogroups, by.x=c("protein_id", "OG_id"), by.y=c("accessions", "OG_id"))
# Map children to parent GO terms
phrogs_long_archaeplastida_go$goId_parent <- ""
for (i in 1:nrow(go_name_to_id_table)) {
  goIds_curr <- c(go_name_to_id_table$goId_parent[i], unlist(strsplit(go_name_to_id_table$goId_children[i], split=",")))
  phrogs_long_archaeplastida_go$goId_parent[which(phrogs_long_archaeplastida_go$goId %in% goIds_curr)] <- paste0(phrogs_long_archaeplastida_go$goId_parent[which(phrogs_long_archaeplastida_go$goId %in% goIds_curr)], ",", go_name_to_id_table$goId_parent[i])
}
phrogs_long_archaeplastida_go <- phrogs_long_archaeplastida_go %>% separate_rows(goId_parent, sep=",")
phrogs_long_archaeplastida_go <- phrogs_long_archaeplastida_go %>% filter(goId_parent != "")
# Remove duplicated protein-GO pairs
phrogs_long_archaeplastida_go <- phrogs_long_archaeplastida_go[!duplicated(cbind(phrogs_long_archaeplastida_go$protein_id, phrogs_long_archaeplastida_go$goId_parent)),]
phrogs_long_archaeplastida_go_counts <- phrogs_long_archaeplastida_go %>% group_by(PROG_id) %>% mutate(n_proteins_total = length(unique(protein_id)), n_GO_proteins = length(unique(protein_id)), n_species = length(unique(taxonId))) %>% group_by(PROG_id, taxonId) %>% mutate(weight = 1 / length(unique(protein_id))) %>% group_by(PROG_id, goId_parent) %>% summarize(n_proteins_total = unique(n_proteins_total), n_GO_proteins = unique(n_GO_proteins), n_species = unique(n_species), sum_weight = sum(weight)) %>% mutate(normalized_weight = sum_weight / n_species)
# Filter by >= 0.5 normalized weight per species to assign GO term
phrogs_long_archaeplastida_go_counts <- phrogs_long_archaeplastida_go_counts %>% filter(normalized_weight >= 0.5)

bl_stem <- bl_df %>% filter(!is.na(neglog10_normalized_stem_bl) & !is.infinite(neglog10_normalized_stem_bl))
bl_dup <- bl_df %>% filter(!is.na(neglog10_normalized_oldest_duplication_bl) & !is.infinite(neglog10_normalized_oldest_duplication_bl))
bl_dup <- bl_dup %>% filter(is.na(neglog10_normalized_stem_bl))
bl_combined <- data.frame(OG_id = c(bl_stem$OG_id, bl_dup$OG_id), parent_PhROG_id = c(bl_stem$parent_PhROG_id, bl_dup$parent_PhROG_id), child_PhROG_id = c(bl_stem$child_PhROG_id, bl_dup$child_PhROG_id), type = c(rep("Acquisition", nrow(bl_stem)), rep("Duplication", nrow(bl_dup))), is_MitoCarta = c(bl_stem$is_MitoCarta, bl_dup$is_MitoCarta), neglog10_normalized_branch_length = c(bl_stem$neglog10_normalized_stem_bl, bl_dup$neglog10_normalized_oldest_duplication_bl))
bl_go <- merge(bl_combined, phrogs_long_archaeplastida_go_counts, by.x="child_PhROG_id", by.y="PROG_id", all.x=TRUE)
# bl_go <- merge(bl_df, phrogs_long_archaeplastida_go_counts, by.x="child_PhROG_id", by.y="PROG_id")
bl_go$GO_name <- go_name_to_id_table$name[match(bl_go$goId_parent, go_name_to_id_table$goId_parent)]
temp <- bl_go %>% filter(!is.na(neglog10_normalized_branch_length) & !is.infinite(neglog10_normalized_branch_length))
temp <- temp %>% filter(!GO_name %in% c("mitochondrion", "cilium", "mitochondrial ribosome", "cell wall", "nuclear envelope"))
bl_MitoCarta <- bl_go %>% filter(is_MitoCarta) %>% filter(!is.na(neglog10_normalized_branch_length) & !is.infinite(neglog10_normalized_branch_length))
bl_MitoCarta$GO_name <- "mitochondrion"
temp <- rbind(temp[,colnames(temp) %in% colnames(bl_MitoCarta)], bl_MitoCarta)
temp <- temp %>% filter(!is.na(GO_name))
temp <- temp %>% group_by(GO_name) %>% filter(!duplicated(child_PhROG_id))
temp_median <- temp %>% group_by(GO_name) %>% summarize(median = median(neglog10_normalized_branch_length))
temp_median <- temp_median[order(temp_median$median, decreasing=TRUE),]
temp$GO_name <- factor(temp$GO_name, levels = c(temp_median$GO_name))
temp_count <- temp %>% group_by(GO_name) %>% summarize(n = n())

# All-vs-all pairwise Wilcoxon rank-sum test, BH correction
pairwise.wilcox.test(temp$neglog10_normalized_branch_length, temp$GO_name, p.adjust.method = "BH", exact = FALSE)

# pdf('combined_length_GO_cellular_compartments_Archaeplastida.pdf', height = 4, width = 4)
# ggplot(temp, aes(x = neglog10_normalized_branch_length, y = GO_name, fill = GO_name)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") + geom_text(
#     data = temp_count,
#     aes(x = Inf, label = n),
#     hjust = 1, # Align text to the right
#     nudge_y = 0.25, # Move the labels slightly above the ridgeline
#     size = 2,
#     color = "black" 
#   ) + xlim(-1.5,3)
# dev.off()






