suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
original_fasta_filename <- args[1]
renamed_fasta_filename <- args[2]

original_fasta <- readAAStringSet(original_fasta_filename)

## Reformat accessions for original fasta files
renamed_fasta <- original_fasta

# Map lta2019 to ltaref
# Read in Lta mapping for Lta2019 to Ltaref
lta_mapping <- read.delim(here("data/orthogroups/idmapping", "map.ltaref.lta.exact.txt"), header=FALSE)
colnames(lta_mapping) <- c("ltaref_id", "lta2019_id")
lta_mapping$ltaref_id <- paste0("5689_", lta_mapping$ltaref_id)
lta_mapping$lta2019_id[lta_mapping$lta2019_id == ""] <- NA
lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)] <- paste0("5689_", lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
lta_mapping <- lta_mapping %>% filter(!is.na(lta2019_id))
names(renamed_fasta)[names(renamed_fasta) %in% lta_mapping$lta2019_id] <- lta_mapping$ltaref_id[match(names(renamed_fasta)[names(renamed_fasta) %in% lta_mapping$lta2019_id], lta_mapping$lta2019_id)]

# Change colons to underscores
names(renamed_fasta) <- gsub(":", "_", names(renamed_fasta), fixed=TRUE)
# Rename incorrectly assigned mtDNA taxids
names(renamed_fasta) <- gsub("^192875_", "595528_", names(renamed_fasta))

writeXStringSet(renamed_fasta, file=renamed_fasta_filename)
