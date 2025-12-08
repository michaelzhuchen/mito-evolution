library(tidyverse)
library(Biostrings)
library(ape)
library(RColorBrewer)
library(basetheme)

### Concatenate individual protein msa files into a supermatrix

data_dir <- here("data/species_phylogeny/single_protein_msa")
msa_list <- read.csv(paste0(data_dir, "/trim_msa_list.txt"), header=FALSE)$V1
# Get list of all species taxids
all_species_taxids <- read.table(here("data/species_phylogeny/guide_tree", "cytosolic_ribosomal_proteins_97.5pct.spp_euk.tophit.from.archaea_euk657spp.prokspp_other.opisthokonta_parasitic.plants_BaSk_CRuMs_bacteria.nomissing.taxlevel6_combined_taxids.txt"), header=FALSE)$V1
exclude_taxids <- c(71139, 112509, 1076696, 158149)
all_species_taxids <- all_species_taxids[!all_species_taxids %in% c(exclude_taxids)] # Remove low quality species with <=12 ribosomal protein homologs
all_species_taxid_names <- all_species_taxids
# Reformat prokgroups names
all_species_taxid_names[all_species_taxid_names %in% uniprot_proteomes_tax$TaxId[uniprot_proteomes_tax$domain %in% c("Bacteria", "Archaea")]] <- uniprot_proteomes_tax$ScientificName[uniprot_proteomes_tax$domain %in% c("Bacteria", "Archaea")]
all_species_taxid_names <- gsub(" ", "", all_species_taxid_names)
all_species_taxid_names <- gsub("\\/", "", all_species_taxid_names)


# Function to count occurrences of a character in a string
count_occurrences <- function(input_string, character_to_count) {
  matches <- gregexpr(character_to_count, input_string)
  return(length(unlist(matches)))
}

missing_species_list <- c()
start_index <- c()
end_index <- c()

for (gene_msa in msa_list) {
  print(gene_msa)
  msa <- readAAStringSet(paste0(data_dir, "/", gene_msa))
  
  # Extract taxid
  names(msa) <- sub(pattern="_.*", replacement="", names(msa)) # Uniprot proteomes
  
  # Identify missing species and add in blank sequences for them
  missing_species <- all_species_taxid_names[which(!(all_species_taxid_names %in% names(msa)))]
  missing_species_list <- c(missing_species_list, missing_species)
  empty_seq <- AAStringSet(paste0(strrep("-", width(msa)[1])))
  for (missing_id in missing_species) {
    names(empty_seq) <- missing_id
    msa <- c(msa, empty_seq)
  }
  
  # Select the top hit from each species based on which paralog has the fewest gaps
  dup_taxids <- unique(names(msa)[which(duplicated(names(msa)))])
  fewest_gaps_list <- c()
  for (dup_taxid in dup_taxids) {
    msa_indices <- which(names(msa) == dup_taxid)
    gap_counts <- sapply(msa[msa_indices], count_occurrences, "-")
    fewest_gaps_for_species <- msa_indices[which(gap_counts == min(gap_counts))[1]] # Break ties by taking the first one
    fewest_gaps_list <- c(fewest_gaps_list, fewest_gaps_for_species)
  }
  # Get a single gene for each species by combining the single copy genes and multicopy top hit sequences
  single_copy_msa <- msa[which(!names(msa) %in% dup_taxids)]
  fewest_gaps_multicopy_msa <- msa[fewest_gaps_list]
  top_hit_msa <- c(single_copy_msa, fewest_gaps_multicopy_msa)
  
  # select species and reorder msa according to the species taxids
  top_hit_msa_reorder <- top_hit_msa[as.character(all_species_taxid_names)]
  
  # concatenate msa
  if (gene_msa == msa_list[1]) {
    concat_msa <- top_hit_msa_reorder
    start_index <- 1
    end_index <- width(top_hit_msa_reorder)[1]
  } else {
    concat_msa <- xscat(concat_msa, top_hit_msa_reorder)
    new_start_index <- tail(end_index, 1) + 1
    start_index <- c(start_index, new_start_index)
    end_index <- c(end_index, width(concat_msa)[1])
  }
}

names(concat_msa) <- all_species_taxids


### Generate NEXUS partition file
gene_ids <- gsub("\\..*", "", msa_list)
nexus_file_lines <- c("#nexus", "begin sets;")
partitions <- paste0("\tcharset ", gene_ids, " = ", start_index, "-", end_index, ";")
nexus_file_lines <- c(nexus_file_lines, partitions)
nexus_file_lines <- c(nexus_file_lines, "end;")

## Write out
# writeXStringSet(concat_msa, paste0(data_dir, "/concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa"))
# writeLines(nexus_file_lines, con = paste0(data_dir, "/", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_partitionfile.nex"))