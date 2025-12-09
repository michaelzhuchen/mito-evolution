suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(castor))

source(here("orthogroups/helpers", "species_overlap_utils.R"))

# Get OG id from command line
args <- commandArgs(trailingOnly=TRUE)
data_dir <- args[1]
outdir <- args[2]
OG_id <- args[3]
suffix <- args[4]
BOOL_STRUCTURE_MERGE <- as.logical(args[5])

## Read in hhsearch results
hhsearch_tsv_file <- paste0(data_dir, "/", OG_id, suffix, ".tsv")

outputfile <- paste0(outdir, "/", OG_id, suffix, "_filter_overlap_coef.tsv")

# Check if file exists and is not empty
if (!file.exists(hhsearch_tsv_file) | file.size(hhsearch_tsv_file) == 0) {
  # Write an empty output file to show completion
  bool_file_create <- file.create(outputfile)
  quit(save = "no")
}

# Analyze hhsearch all vs all OG results for merging splintered OGs
hhsearch_OG <- read.table(hhsearch_tsv_file, sep="\t", header=FALSE)

hhsearch_OG$query_OG <- gsub("\\.fa.*", "", hhsearch_OG$V1)
hhsearch_OG$target_OG <- gsub("\\.fa.*", "", hhsearch_OG$V3)
hhsearch_OG$expect <- hhsearch_OG$V5
hhsearch_OG$score <- hhsearch_OG$V7
hhsearch_OG$prob <- hhsearch_OG$V4

if (BOOL_STRUCTURE_MERGE) {
  hhsearch_OG_filter <- hhsearch_OG %>% filter(query_OG != target_OG)
} else {
  hhsearch_OG_filter <- hhsearch_OG %>% filter(query_OG != target_OG) %>% filter(expect < 1e-3)
}

# Check if table is empty
if (nrow(hhsearch_OG_filter) == 0) {
  # Write an empty output file to show completion
  bool_file_create <- file.create(outputfile)
  quit(save = "no")
}

hhsearch_OG_filter <- hhsearch_OG_filter[order(hhsearch_OG_filter$expect, 1/hhsearch_OG_filter$score, decreasing=FALSE),]
hhsearch_OG_filter <- hhsearch_OG_filter[!duplicated(hhsearch_OG_filter[,c("query_OG", "target_OG")]),]

## Create copies matrix
if (BOOL_STRUCTURE_MERGE) {
  # For hmm+structure merge, directly use the copies matrix after the initial merge
  copies_mat <- readRDS(here("data/orthogroups/first_merge", "OG_all_merged_copies_mat.rds"))
} else {
  # For initial hmm-hmm merge, reformat the original copies table into a copies matrix
  orthogroups_copies <- readRDS(here("data/orthogroups/raw_orthogroups/Orthogroups", "Orthogroups.GeneCount_ACABI.best.isoform.rds"))
  copies_mat <- orthogroups_copies[,2:ncol(orthogroups_copies)]
  rownames(copies_mat) <- orthogroups_copies$Orthogroup
}

## Compute overlap coefficient for all hits
copies_mat_update <- copies_mat
overlap_coef_list <- c()
for (i in 1:nrow(hhsearch_OG_filter)) {
  curr_hhsearch <- hhsearch_OG_filter[i,]
  
  # Grep to find the matching rows in copies matrix, taking word boundaries into account (in case of merged OGs that are comma-separated)
  query_OG_index <- grep(paste0("\\b", curr_hhsearch$query_OG, "\\b"), rownames(copies_mat_update))
  target_OG_index <- grep(paste0("\\b", curr_hhsearch$target_OG, "\\b"), rownames(copies_mat_update))
  
  overlap_coef <- fraction_overlap(copies_mat_update[query_OG_index,], copies_mat_update[target_OG_index,])
  
  overlap_coef_list <- c(overlap_coef_list, overlap_coef)
}

# Combine with hhsearch table and select specific columns to reduce data size
if (BOOL_STRUCTURE_MERGE) {
  hhsearch_OG_filter_overlap_coef <- cbind(hhsearch_OG_filter[,c("query_OG", "target_OG", "expect", "score")], overlap_coef_list, hhsearch_OG_filter[,c("prob")])
} else {
  hhsearch_OG_filter_overlap_coef <- cbind(hhsearch_OG_filter[,c("query_OG", "target_OG", "expect", "score")], overlap_coef_list)
}

## Write out
# write.table(hhsearch_OG_filter_overlap_coef, outputfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

