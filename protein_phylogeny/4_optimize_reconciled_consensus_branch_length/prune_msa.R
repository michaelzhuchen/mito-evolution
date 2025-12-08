# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

# Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
msa_base_dir <- args[1]
origin_base_dir <- args[2]
out_dir <- args[3]
OG_id <- args[4]
BOOL_VERBOSE <- FALSE

dir.create(out_dir, showWarnings = FALSE)

msa_filename <- paste0(msa_base_dir, "/", OG_id, ".faa_clipkit.gappy.msa")
msa <- readAAStringSet(msa_filename)

origin_table <- read.table(here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t")
colnames(origin_table) <- c("OG_id", "n_euk_species_largest_euk_clade", "n_prok_species_largest_prok_clade", "LCA_node", "LCA_node_euks", "origin_domain", "dropped_tips")
origin_table <- origin_table %>% filter(OG_id == OG_id)

if (!is.na(origin_table$dropped_tips)) {
  dropped_tips_list <- unlist(strsplit(origin_table$dropped_tips, split=","))
  if (BOOL_VERBOSE) {
    print(dropped_tips_list)
  }
  msa_pruned <- msa[which(!names(msa) %in% dropped_tips_list),]
} else {
  if (BOOL_VERBOSE) {
    print("No dropped tips")
  }
  msa_pruned <- msa
}

msa_pruned_filename <- paste0(out_dir, "/", OG_id, ".faa_clipkit.gappy_pruned.msa")

writeXStringSet(msa_pruned, msa_pruned_filename)
