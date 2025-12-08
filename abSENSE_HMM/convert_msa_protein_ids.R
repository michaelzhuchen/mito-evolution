suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
original_msa_filename <- args[1]
renamed_msa_filename <- args[2]

original_msa <- readAAStringSet(original_msa_filename)

# Map lta2019 to ltaref
# Read in Lta mapping for Lta2019 to Ltaref
lta_mapping <- read.delim(here("data/orthogroups/idmapping", "map.ltaref.lta.exact.txt"), header=FALSE)
colnames(lta_mapping) <- c("ltaref_id", "lta2019_id")
lta_mapping$ltaref_id <- paste0("5689_", lta_mapping$ltaref_id)
lta_mapping$lta2019_id[lta_mapping$lta2019_id == ""] <- NA
lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)] <- paste0("5689_", lta_mapping$lta2019_id[!is.na(lta_mapping$lta2019_id)])
lta_mapping <- lta_mapping %>% filter(!is.na(lta2019_id))

renamed_msa <- original_msa
names(renamed_msa)[names(renamed_msa) %in% lta_mapping$lta2019_id] <- lta_mapping$ltaref_id[match(names(renamed_msa)[names(renamed_msa) %in% lta_mapping$lta2019_id], lta_mapping$lta2019_id)]

writeXStringSet(renamed_msa, file=renamed_msa_filename)
