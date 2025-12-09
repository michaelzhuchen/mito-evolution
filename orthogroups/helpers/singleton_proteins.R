### Get singleton proteins
tophit_singleton <- read.table(here("data/orthogroups/hmmsearch_singleton_proteins", "hmmsearch_singleton_proteins_expect1e-3_DBSIZE.379668_tophit.tsv"))
colnames(tophit_singleton) <- c("protein_id", "OG_id")

# Split ties into separate rows
tophit_singleton_separaterows <- tophit_singleton %>% separate_rows(OG_id, sep=",")

# Remove obsolete GiardiaDB44 proteins
tophit_singleton_separaterows <- tophit_singleton_separaterows %>% filter(!grepl("^5741_GL50803_", protein_id))

## Remap accessions for ACANB, Tbr, Lta.
# Keep all ACABI transcript IDs that are identical in ACANB and are the best isoform
ACABI_singleton_protein_ids <- tophit_singleton_separaterows$protein_id[grep("^1257118_", tophit_singleton_separaterows$protein_id)]
remove_ACABI_singleton_protein_ids <- ACABI_singleton_protein_ids[which(!ACABI_singleton_protein_ids %in% ACABI_ids_identical_to_ACANB_best_transcript_isoforms_accessions)]
tophit_singleton_separaterows <- tophit_singleton_separaterows[which(!tophit_singleton_separaterows$protein_id %in% remove_ACABI_singleton_protein_ids),]
# Replace ACABI with ACANB ids
tophit_singleton_separaterows$protein_id[which(tophit_singleton_separaterows$protein_id %in% ACANB_ACABI_transcripts$ACABI_transcript_id)] <- ACANB_ACABI_transcripts$ACANB_transcript_id[match(tophit_singleton_separaterows$protein_id[which(tophit_singleton_separaterows$protein_id %in% ACANB_ACABI_transcripts$ACABI_transcript_id)], ACANB_ACABI_transcripts$ACABI_transcript_id)]

# Keep Tbr transcript ids that map to best isoform
tbr_singleton_protein_ids <- tophit_singleton_separaterows$protein_id[grep("^185431_", tophit_singleton_separaterows$protein_id)]
remove_tbr_singleton_protein_ids <- tbr_singleton_protein_ids[which(!tbr_singleton_protein_ids %in% tbr_mapping$transcriptID)]
tophit_singleton_separaterows <- tophit_singleton_separaterows[which(!tophit_singleton_separaterows$protein_id %in% remove_tbr_singleton_protein_ids),]
# Replace Tbr transcript id with best isoform gene id
tophit_singleton_separaterows$protein_id[which(tophit_singleton_separaterows$protein_id %in% tbr_mapping$transcriptID)] <- tbr_mapping$ID[match(tophit_singleton_separaterows$protein_id[which(tophit_singleton_separaterows$protein_id %in% tbr_mapping$transcriptID)], tbr_mapping$transcriptID)]

# Keep lta2019 proteins ids that map to ltaref
lta_singleton_protein_ids <- tophit_singleton_separaterows$protein_id[grep("^5689_", tophit_singleton_separaterows$protein_id)]
remove_lta_singleton_protein_ids <- lta_singleton_protein_ids[which(!lta_singleton_protein_ids %in% lta_mapping$lta2019_id)]
tophit_singleton_separaterows <- tophit_singleton_separaterows[which(!tophit_singleton_separaterows$protein_id %in% remove_lta_singleton_protein_ids),]
# Replace lta2019 protein id with ltaref
tophit_singleton_separaterows$protein_id[which(tophit_singleton_separaterows$protein_id %in% lta_mapping$lta2019_id)] <- lta_mapping$ltaref_id[match(tophit_singleton_separaterows$protein_id[which(tophit_singleton_separaterows$protein_id %in% lta_mapping$lta2019_id)], lta_mapping$lta2019_id)]

# Get summary of new proteins by OG to add to their new OGs in block below
tophit_singleton_summary <- tophit_singleton_separaterows %>% group_by(OG_id) %>% summarize(accessions = paste(unique(protein_id), sep="", collapse=", "))
