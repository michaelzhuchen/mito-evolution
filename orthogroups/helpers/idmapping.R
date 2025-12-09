### Get ID mapping tables for updated protein ids

# Read in ACABI->ACANB mapping
ACABI_best_transcript_isoforms <- read.table(here("data/orthogroups/idmapping", "ACABI_best_transcript_id_list.txt"))$V1
ACANB_best_transcript_isoforms <- read.delim(here("data/orthogroups/idmapping", "ACANBgene.genes_best.isoform.txt"))$transcriptID
ACABI_all_transcript_isoforms <- read.table(here("data/orthogroups/idmapping", "ACABI_all_transcript_id_list.txt"))$V1
ACANB_ACABI_mapping <- read.csv(here("data/orthogroups/idmapping", "map.ACANB.ACABI.ACANBI.csv"))
ACABI_ids_identical_to_ACANB_best_transcript_isoforms_accessions <- paste0("1257118_", ACANB_ACABI_mapping$ACABI_transcriptID[ACANB_ACABI_mapping$ACANB_vs_ACABI_same_protein == 1 & ACANB_ACABI_mapping$ACANB_transcriptID %in% ACANB_best_transcript_isoforms])
ACANB_ACABI_transcripts <- data.frame(ACABI_transcript_id = ACANB_ACABI_mapping$ACABI_transcriptID, ACANB_transcript_id = ACANB_ACABI_mapping$ACANB_transcriptID)
ACANB_ACABI_transcripts <- ACANB_ACABI_transcripts %>% rowwise() %>% mutate(ACABI_transcript_id = paste0("1257118_", ACABI_transcript_id), ACANB_transcript_id = paste0("1257118_", ACANB_transcript_id))

# Read in new Tbr mapping for fasta with best isoform per gene
tbr_mapping <- read.table(here("data/orthogroups/idmapping", "tbrgene.mapping.txt"), header=TRUE)
tbr_mapping$transcriptID <- gsub(":", "_", tbr_mapping$transcriptID, fixed=TRUE)
tbr_mapping$ID <- paste0("185431_", tbr_mapping$ID)
tbr_mapping$transcriptID <- paste0("185431_", tbr_mapping$transcriptID)

# Read in new Lta mapping for Lta2019 to Ltaref
lta_mapping <- read.delim(here("data/orthogroups/idmapping", "map.ltaref.lta.exact.txt"), header=FALSE)
colnames(lta_mapping) <- c("ltaref_id", "lta2019_id")
lta_mapping$ltaref_id <- paste0("5689_", lta_mapping$ltaref_id)
lta_mapping$lta2019_id[lta_mapping$lta2019_id == ""] <- NA
lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)] <- paste0("5689_", lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
