### Get initially merged orthogroups (non-singletons)

# Load libraries
library(here)
library(tidyverse)
library(ape)
library(castor)
library(Biostrings)

### Read in mapping for updated IDs
source(here("orthogroups/helpers", "idmapping.R"))

### Build initially merged orthogroups

# Read in raw orthogroups
orthogroups <- read.delim(here("data/orthogroups/raw_orthogroups/Orthogroups", "Orthogroups.tsv"), header=TRUE)

# Read in hmm-merged OGs
orthogroups_merge <- read.table(here("data/orthogroups/first_merge", "first_merged_orthogroups.tsv"), sep="\t", header = TRUE)
# Assign new OG ids for merged OGs
new_OG_ids <- sprintf("MOG%07d", 1:nrow(orthogroups_merge))
orthogroups_merge$Orthogroup_original <- orthogroups_merge$Orthogroup
orthogroups_merge$Orthogroup[grep(",", orthogroups_merge$Orthogroup)] <- new_OG_ids[grep(",", orthogroups_merge$Orthogroup)]

# Add merged OGs to the raw OGs, replacing raw OGs with merged OGs
all_merged_OG_ids <- unlist(strsplit(orthogroups_merge$Orthogroup_original, split=", "))
orthogroups <- orthogroups[which(!orthogroups$Orthogroup %in% all_merged_OG_ids),]
orthogroups$Orthogroup_original <- orthogroups$Orthogroup
orthogroups <- rbind(orthogroups_merge, orthogroups)

# Remove excluded taxa that have <=8 cytoribo proteins or obsolete annotation version (GiardiaDB44)
exclude_taxids <- c(71139, 112509, 1076696, 158149, 5741)
exclude_taxids_colnames <- paste0("X", exclude_taxids)
orthogroups <- orthogroups[,which(!colnames(orthogroups) %in% exclude_taxids_colnames)]

### Map old accessions to new accessions for Aca, Tbr, Lta.
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

# Collapse the OG ID and protein IDs for sequence retrieval
orthogroups$collapse_accessions <- apply(orthogroups[,colnames(orthogroups) != "Orthogroup_original"], 1, function(row) paste(str_trim(row[row != ""]), sep="", collapse = ";"))
orthogroups$collapse_accessions <- gsub(";", " >", orthogroups$collapse_accessions)
orthogroups$collapse_accessions <- gsub(", ", " >", orthogroups$collapse_accessions)
orthogroups$collapse_accessions <- gsub(",", " >", orthogroups$collapse_accessions)

## Write out
# write.table(orthogroups[,c("Orthogroup", "Orthogroup_original")], here("data/orthogroups/first_merge", "initially_merged_OG_id_mapping.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

