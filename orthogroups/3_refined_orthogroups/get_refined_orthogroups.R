### Get refined orthogroups

# Load libraries
library(here)
library(tidyverse)
library(ape)
library(castor)
library(Biostrings)

## Start from raw orthogroups
# Read in raw orthogroups
orthogroups <- read.delim(here("data/orthogroups/raw_orthogroups/Orthogroups", "Orthogroups.tsv"), header=TRUE)

## First and second merges
# Read in hmm-merged OGs
orthogroups_merge <- read.table(here("data/orthogroups/first_merge", "first_merged_orthogroups.tsv"), sep="\t", header = TRUE)
# Assign new OG ids for merged OGs
new_OG_ids <- sprintf("MOG%07d", 1:nrow(orthogroups_merge))
orthogroups_merge$Orthogroup_original <- orthogroups_merge$Orthogroup
orthogroups_merge$Orthogroup[grep(",", orthogroups_merge$Orthogroup)] <- new_OG_ids[grep(",", orthogroups_merge$Orthogroup)]

# Read in hmm+structure merged OGs
orthogroups_foldseek_merge <- read.table(here("data/orthogroups/second_merge", "second_merged_orthogroups.tsv"), sep="\t", header = TRUE)
# Assign new OG ids for hmm+structure merged OGs
new_OG_ids <- sprintf("SOG%07d", 1:nrow(orthogroups_foldseek_merge))
orthogroups_foldseek_merge$Orthogroup_original <- orthogroups_foldseek_merge$Orthogroup
orthogroups_foldseek_merge$Orthogroup[grep(",", orthogroups_foldseek_merge$Orthogroup)] <- new_OG_ids[grep(",", orthogroups_foldseek_merge$Orthogroup)]

# Remove redundant hmm-merged OGs that are further merged by hmm+structure
copies_mat_update <- readRDS(here("data/orthogroups/second_merge", "second_merged_copies_mat.rds"))
copies_mat_update <- copies_mat_update[!grepl("X", rownames(copies_mat_update)),]
merged_OG_ids <- rownames(copies_mat_update)[grep(",", rownames(copies_mat_update))]
merged_OG_df <- data.frame(OG_index=1:length(merged_OG_ids), OG_id=merged_OG_ids)
merged_OG_seprows <- merged_OG_df %>% separate_rows(OG_id, sep=",")
orthogroups_merge <- orthogroups_merge[which(!orthogroups_merge$Orthogroup %in% merged_OG_seprows$OG_id),]

# Add merged OGs to the raw OGs, replacing raw OGs with merged OGs
all_merged_OG_ids <- unlist(strsplit(orthogroups_merge$Orthogroup_original, split=", "))
all_foldseek_merged_OG_ids <- unlist(strsplit(orthogroups_foldseek_merge$Orthogroup_original, split=", "))
orthogroups <- orthogroups[which(!orthogroups$Orthogroup %in% c(all_merged_OG_ids, all_foldseek_merged_OG_ids)),]
orthogroups$Orthogroup_original <- orthogroups$Orthogroup
orthogroups <- rbind(orthogroups_merge, orthogroups)
orthogroups <- rbind(orthogroups_foldseek_merge, orthogroups)

# Remove excluded taxa that have <=8 cytoribo proteins or obsolete annotation version (GiardiaDB44)
exclude_taxids <- c(71139, 112509, 1076696, 158149, 5741)
exclude_taxids_colnames <- paste0("X", exclude_taxids)
orthogroups <- orthogroups[,which(!colnames(orthogroups) %in% exclude_taxids_colnames)]


### Map old accessions to new accessions for Aca, Tbr, Lta.
## Read in mapping for updated IDs
source(here("orthogroups/helpers", "idmapping.R"))

## Map ACABI->ACANB and keep best transcript isoforms.
# Keep all ACABI transcript IDs that are identical in ACANB and are the best isoform
filter_items <- function(item_string, allowed_items) {
  items <- unlist(strsplit(item_string, ", "))
  filtered_items <- intersect(items, allowed_items)
  paste(filtered_items, collapse = ", ")
}
orthogroups$X1257118 <- sapply(orthogroups$X1257118, filter_items, ACABI_ids_identical_to_ACANB_best_transcript_isoforms_accessions)

# Replace ACABI transcript ids with ACANB ids
replace_items <- function(item_string, mapping_table) {
  items <- unlist(strsplit(item_string, ", "))
  mapped_items <- mapping_table$ACANB_transcript_id[which(mapping_table$ACABI_transcript_id %in% items)]
  paste(mapped_items, collapse = ", ")
}
orthogroups$X1257118 <- sapply(orthogroups$X1257118, replace_items, ACANB_ACABI_transcripts)

## Map new Tbr accessions and keep only the best transcript isoforms for Tbr (185431)
# Keep all Tbr transcript IDs that are identical in new annotation
orthogroups$X185431 <- sapply(orthogroups$X185431, filter_items, tbr_mapping$transcriptID)
# Replace Tbr transcript IDs with the best isoform gene id
replace_items <- function(item_string, mapping_table) {
  items <- unlist(strsplit(item_string, ", "))
  mapped_items <- mapping_table$ID[which(mapping_table$transcriptID %in% items)]
  paste(mapped_items, collapse = ", ")
}
orthogroups$X185431 <- sapply(orthogroups$X185431, replace_items, tbr_mapping)

## Map new Lta accessions and keep only the identical proteins in ltaref (5689)
# Keep all Lta proteins that are identical in ltaref
orthogroups$X5689 <- sapply(orthogroups$X5689, filter_items, lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
# Replace Lta2019 protein ids with ltaref protein ids
replace_items <- function(item_string, mapping_table) {
  items <- unlist(strsplit(item_string, ", "))
  mapped_items <- mapping_table$ltaref_id[which(mapping_table$lta2019_id %in% items)]
  paste(mapped_items, collapse = ", ")
}
orthogroups$X5689 <- sapply(orthogroups$X5689, replace_items, lta_mapping)

colnames(orthogroups) <- gsub("^X", "", colnames(orthogroups))


## Add in singleton genes
source(here("orthogroups/helpers", "singleton_proteins.R"))
# orthogroups <- add_proteins_to_orthogroups(orthogroups, tophit_singleton_separaterows)
orthogroups$singleton_protein_accessions <- ""
tophit_singleton_summary <- tophit_singleton_summary[tophit_singleton_summary$OG_id %in% orthogroups$Orthogroup,]
orthogroups$singleton_protein_accessions[match(tophit_singleton_summary$OG_id, orthogroups$Orthogroup)] <- tophit_singleton_summary$accessions

## Add in added species proteins
source(here("orthogroups/helpers", "added_proteins.R"))
# orthogroups <- add_proteins_to_orthogroups(orthogroups, tophit_add_species_separaterows)
orthogroups$added_species_protein_accessions <- ""
tophit_add_species_summary <- tophit_add_species_summary[tophit_add_species_summary$OG_id %in% orthogroups$Orthogroup,]
orthogroups$added_species_protein_accessions[match(tophit_add_species_summary$OG_id, orthogroups$Orthogroup)] <- tophit_add_species_summary$accessions

## Add in fusion genes
source(here("orthogroups/helpers", "fusion_proteins.R"))
# colnames(combined_fusion_protein_filter_disjoint_separaterows_newhit)[2] <- "OG_id"
# combined_fusion_protein_filter_disjoint_separaterows_newhit <- combined_fusion_protein_filter_disjoint_separaterows_newhit[,c("protein_id", "OG_id")]
# orthogroups <- add_proteins_to_orthogroups(orthogroups, combined_fusion_protein_filter_disjoint_separaterows_newhit)
orthogroups$fusion_protein_accessions <- ""
fusion_protein_summary <- fusion_protein_summary[fusion_protein_summary$disjoint_hits %in% orthogroups$Orthogroup,]
orthogroups$fusion_protein_accessions[match(fusion_protein_summary$disjoint_hits, orthogroups$Orthogroup)] <- fusion_protein_summary$accessions

# Collapse the OG ID and protein IDs for sequence retrieval
orthogroups$collapse_accessions <- apply(orthogroups[,colnames(orthogroups) != "Orthogroup_original"], 1, function(row) paste(str_trim(row[row != ""]), sep="", collapse = ";"))
orthogroups$collapse_accessions <- gsub(";", " >", orthogroups$collapse_accessions)
orthogroups$collapse_accessions <- gsub(", ", " >", orthogroups$collapse_accessions)
orthogroups$collapse_accessions <- gsub(",", " >", orthogroups$collapse_accessions)

## Assign OG IDs to singletons
# Generate UOG ids for initial set of singleton proteins from raw OGs
singleton_proteins_raw <- read.table(here("data/orthogroups/raw_orthogroups/Orthogroups", "raw_orthogroup_singleton_protein_accessions.txt"), sep="\t", header=FALSE)
colnames(singleton_proteins_raw) <- "accession"
# Remove obsolete GiardiaDB44 proteins
singleton_proteins_raw <- singleton_proteins_raw %>% filter(!grepl("^5741_GL50803_", accession))
# Assign OG ids to singleton proteins
singleton_proteins_raw$OG_id <- sprintf("UOG%07d", 1:nrow(singleton_proteins_raw))

## Map to updated ids for ACANB, Tbr, and Lta
# Keep all ACABI transcript IDs that are identical in ACANB and are the best isoform
singleton_proteins <- singleton_proteins_raw[,c("OG_id", "accession")]
ACABI_singleton_protein_ids <- singleton_proteins$accession[grep("^1257118_", singleton_proteins$accession)]
remove_ACABI_singleton_protein_ids <- ACABI_singleton_protein_ids[which(!ACABI_singleton_protein_ids %in% ACABI_ids_identical_to_ACANB_best_transcript_isoforms_accessions)]
singleton_proteins <- singleton_proteins[which(!singleton_proteins$accession %in% remove_ACABI_singleton_protein_ids),]
# Replace ACABI with ACANB protein
singleton_proteins$accession[which(singleton_proteins$accession %in% ACANB_ACABI_transcripts$ACABI_transcript_id)] <- ACANB_ACABI_transcripts$ACANB_transcript_id[match(singleton_proteins$accession[which(singleton_proteins$accession %in% ACANB_ACABI_transcripts$ACABI_transcript_id)], ACANB_ACABI_transcripts$ACABI_transcript_id)]

# Keep all Tbr transcript ids that are best isoform in new Tbr annotation
Tbr_singleton_protein_ids <- singleton_proteins$accession[grep("^185431_", singleton_proteins$accession)]
remove_Tbr_singleton_protein_ids <- Tbr_singleton_protein_ids[which(!Tbr_singleton_protein_ids %in% tbr_mapping$transcriptID)]
singleton_proteins <- singleton_proteins[which(!singleton_proteins$accession %in% remove_Tbr_singleton_protein_ids),]
# Replace Tbr transcript id with gene id
singleton_proteins$accession[which(singleton_proteins$accession %in% tbr_mapping$transcriptID)] <- tbr_mapping$ID[match(singleton_proteins$accession[which(singleton_proteins$accession %in% tbr_mapping$transcriptID)], tbr_mapping$transcriptID)]

# Keep lta2019 proteins ids that map to ltaref
lta_singleton_protein_ids <- singleton_proteins$accession[grep("^5689_", singleton_proteins$accession)]
remove_lta_singleton_protein_ids <- lta_singleton_protein_ids[which(!lta_singleton_protein_ids %in% lta_mapping$lta2019_id)]
singleton_proteins <- singleton_proteins[which(!singleton_proteins$accession %in% remove_lta_singleton_protein_ids),]
# Replace lta2019 protein id with ltaref
singleton_proteins$accession[which(singleton_proteins$accession %in% lta_mapping$lta2019_id)] <- lta_mapping$ltaref_id[match(singleton_proteins$accession[which(singleton_proteins$accession %in% lta_mapping$lta2019_id)], lta_mapping$lta2019_id)]

# Get list of singleton proteins with no top hit
singleton_proteins_unmerged <- singleton_proteins[!singleton_proteins$accession %in% tophit_singleton_separaterows$protein_id,]

# Get list of added species proteins with no top hit
protein_ids_81824 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "81824.fasta")))
protein_ids_595528 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "595528.fasta")))
protein_ids_667725 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "667725.fasta")))
protein_ids_691883 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "691883.fasta")))
protein_ids_946362 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "946362.fasta")))
protein_ids_5741 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "5741.fasta")))
protein_ids_1257118 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "ACANB_no.asterisk_new.proteins_rename.pep")))
protein_ids_BaSk <- names(readAAStringSet(here("data/orthogroups/added_fastas", "BaSk_rename_combined.fasta")))
protein_ids_CRuMs <- names(readAAStringSet(here("data/orthogroups/added_fastas", "CRuMs_rename_combined.fasta")))
protein_ids_parasiticplants <- names(readAAStringSet(here("data/orthogroups/added_fastas", "Parasitic_plants_rename_combined.fasta")))
protein_ids_ltaref <- names(readAAStringSet(here("data/orthogroups/added_fastas", "5689_no.asterisk_new.proteins.fasta")))
protein_ids_mtDNA_new_proteins <- names(readAAStringSet(here("data/orthogroups/added_fastas", "mtDNA_protein_fasta_rename_combined.faa")))
protein_ids_mtDNA_new_proteins_2025.04.07 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "mtDNA_protein_fasta_rename_added.2025.04.07_combined.faa")))
protein_ids_mtDNA_new_proteins_2025.05.31 <- names(readAAStringSet(here("data/orthogroups/added_fastas", "mtDNA_protein_fasta_rename_added.2025.05.31_combined.faa")))

# Remove incorrectly matched accessions from mtDNA and rename incorrect taxid
protein_ids_mtDNA_new_proteins <- protein_ids_mtDNA_new_proteins[!grepl("^2086695_", protein_ids_mtDNA_new_proteins)]
protein_ids_mtDNA_new_proteins[grep("^192875_", protein_ids_mtDNA_new_proteins)] <- gsub("^192875_", "595528_", protein_ids_mtDNA_new_proteins[grep("^192875_", protein_ids_mtDNA_new_proteins)])

all_protein_ids_from_added_fastas <- c(protein_ids_81824, protein_ids_595528, protein_ids_667725, protein_ids_691883, protein_ids_946362, protein_ids_5741, protein_ids_1257118, protein_ids_BaSk, protein_ids_CRuMs, protein_ids_parasiticplants, protein_ids_mtDNA_new_proteins, protein_ids_mtDNA_new_proteins_2025.04.07, protein_ids_ltaref, protein_ids_mtDNA_new_proteins_2025.05.31)
all_protein_ids_from_added_fastas <- gsub(":", "_", all_protein_ids_from_added_fastas, fixed=TRUE)

# Get list of added species proteins with no top hit
singleton_protein_ids_from_added_fastas <- all_protein_ids_from_added_fastas[which(!all_protein_ids_from_added_fastas %in% tophit_add_species_separaterows$protein_id)]
singleton_protein_ids_from_added_fastas_df <- data.frame(OG_id = sprintf("UOG%07d", nrow(singleton_proteins_raw) + 1:length(singleton_protein_ids_from_added_fastas)), accession = singleton_protein_ids_from_added_fastas)

singleton_proteins_unmerged <- rbind(singleton_proteins_unmerged, singleton_protein_ids_from_added_fastas_df)

singleton_proteins_unmerged$accession <- paste0(">", singleton_proteins_unmerged$accession)
singleton_proteins_unmerged$collapse_accession <- paste0(singleton_proteins_unmerged$OG_id, " ", singleton_proteins_unmerged$accession)

# Update lists: Join with non-singleton OG ids and protein accessions
orthogroups_collapsed <- data.frame(accession = c(orthogroups$collapse_accessions, singleton_proteins_unmerged$collapse_accession))
orthogroups_collapsed$Orthogroup <- gsub(" .*", "", orthogroups_collapsed$accession)

### Generate long format output. Calculate copy counts per species per OG
ogs_long <- orthogroups_collapsed %>% separate_rows(accession, sep=" >")
ogs_long <- ogs_long[which(ogs_long$accession != ogs_long$Orthogroup),]
ogs_long$taxid <- gsub("_.*", "", ogs_long$accession)

# Mark primary OG
ogs_long_merge <- ogs_long
ogs_long_merge$rowindex <- 1:nrow(ogs_long_merge)
ogs_long_merge <- merge(ogs_long_merge, combined_fusion_protein_filter_disjoint_separaterows_newhit, by.x=c("accession", "Orthogroup"), by.y=c("protein_id", "disjoint_hits"))
ogs_long$BOOL_PRIMARY_OG <- TRUE
ogs_long$BOOL_PRIMARY_OG[ogs_long_merge$rowindex] <- FALSE

### Generate wide format output
ogs_wide <- ogs_long %>% select(-BOOL_PRIMARY_OG) %>% pivot_wider(names_from = taxid, values_from = accession, values_fn = list(accession = ~paste(., collapse = ",")), values_fill = "")

# Reorder by taxonomy
uniprot_proteomes_673euks_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
uniprot_proteomes_673euks_tax <- uniprot_proteomes_673euks_tax[order(uniprot_proteomes_673euks_tax$Lineage),]
ogs_wide <- cbind(ogs_wide$Orthogroup, ogs_wide[,match(uniprot_proteomes_673euks_tax$tree_id, colnames(ogs_wide))])

# Rename to species names
colnames(ogs_wide) <- uniprot_proteomes_673euks_tax$ScientificName[match(colnames(ogs_wide), uniprot_proteomes_673euks_tax$tree_id)]
colnames(ogs_wide)[1] <- "Orthogroup_ID"

### Write out
# write.table(ogs_long, here("data/orthogroups/refined_orthogroups", "refined_OGs_euk673spp_long.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
# write.table(ogs_wide, here("data/orthogroups/refined_orthogroups", "refined_OGs_euk673spp_wide.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
