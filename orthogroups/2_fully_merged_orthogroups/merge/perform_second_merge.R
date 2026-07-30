### Perform second merge based on hmm-hmm alignment, structure alignment, and species overlap

# Load libraries
library(ape)
library(tidyverse)
library(castor)

source(here("orthogroups/helpers", "species_overlap_utils.R"))

suffix <- "second_merge"
outdir <- here("data/orthogroups/second_merge")

print("Reading in data")

## Read in hhsearch data
hhsearch_OG <- read.table(here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_overlap_coef0.3_bidir_top50hits_mean.bitscore.prob_combined_processed.tsv"), sep="\t")
colnames(hhsearch_OG) <- c("query_OG", "target_OG", "prob", "score", "overlap_coef")

hhsearch_OG_filter <- hhsearch_OG %>% filter(query_OG != target_OG) %>% filter(prob > 0.9) # for hmm+structure merge
hhsearch_OG_filter <- hhsearch_OG_filter[order(hhsearch_OG_filter$prob, hhsearch_OG_filter$score, decreasing=TRUE),]
hhsearch_OG_filter <- hhsearch_OG_filter[!duplicated(hhsearch_OG_filter[,c("query_OG", "target_OG")]),]

# Get top hit per query for nonbidir merge
hhsearch_OG_filter_bidir <- hhsearch_OG_filter

# Read in the full updated copies matrix after hmm-hmm merge 
copies_mat <- readRDS(here("data/orthogroups/first_merge", "OG_all_merged_copies_mat.rds"))

# Get subset of OGs that will potentially be merged
potential_merge_OGs <- unique(c(hhsearch_OG_filter_bidir$query_OG, hhsearch_OG_filter_bidir$target_OG, hhsearch_OG_filter_nonbidir$query_OG, hhsearch_OG_filter_nonbidir$target_OG))
copies_mat_update <- copies_mat[which(rownames(copies_mat) %in% potential_merge_OGs),]

print(paste0("Total of ", nrow(copies_mat_update), " unique OGs with candidate merges..."))

## Merge OGs with their bidir and nonbidir hits in a greedy way, taking species overlap coefficient into account
# Do this for bidirectional hits first.
print(paste0("Processing ", nrow(hhsearch_OG_filter_bidir), " bidir hits..."))
fraction_overlap_threshold <- 0.3

for (i in 1:nrow(hhsearch_OG_filter_bidir)) {
  curr_hhsearch <- hhsearch_OG_filter_bidir[i,]

  # Grep to find the matching rows in copies matrix, taking word boundaries into account (for merged OGs that are comma-separated)
  query_OG_index <- grep(paste0("\\b", curr_hhsearch$query_OG, "\\b"), rownames(copies_mat_update))
  target_OG_index <- grep(paste0("\\b", curr_hhsearch$target_OG, "\\b"), rownames(copies_mat_update))

  # If already merged
  if (query_OG_index == target_OG_index) {
    next
  }

  overlap_coef <- fraction_overlap(copies_mat_update[query_OG_index,], copies_mat_update[target_OG_index,])

  # Merge the OGs
  if (overlap_coef < fraction_overlap_threshold) {
    copies_mat_update[target_OG_index,] <- colSums(copies_mat_update[c(query_OG_index,target_OG_index),])
    query_OG_name <- rownames(copies_mat_update)[query_OG_index]
    target_OG_name <- rownames(copies_mat_update)[target_OG_index]
    rownames(copies_mat_update)[target_OG_index] <- paste0(target_OG_name, ",", query_OG_name)

    # copies_mat_update <- copies_mat_update[-query_OG_index,]
    rownames(copies_mat_update)[query_OG_index] <- paste0("X", query_OG_index)

    print(paste0("Bidir merge (", i, " of ", nrow(hhsearch_OG_filter_bidir), ")"))

  }
}

# saveRDS(copies_mat_update, here("data/orthogroups/second_merge", "second_merged_copies_mat.rds"))

