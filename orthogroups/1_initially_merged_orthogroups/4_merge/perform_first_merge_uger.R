### Perform first merge based on hmm-hmm alignment and species overlap

# Load libraries
library(ape)
library(tidyverse)
library(castor)

source(here("orthogroups/helpers", "species_overlap_utils.R"))

suffix <- "first_merge"
outdir <- here("data/orthogroups/first_merge")

print("Reading in data")

## Preprocess OG copy counts
# Read from RDS
orthogroups_copies <- readRDS(here("data/orthogroups/raw_orthogroups/Orthogroups", "Orthogroups.GeneCount_ACABI.best.isoform.rds"))

## Read in hhsearch data
# Read in hhsearch + species overlap coef data, identify bidirectional hits and nonbidirectional hits
hhsearch_OG <- read.table(here("data/orthogroups/first_merge", "combined_OG_all_vs_all_hhsuite_expect1e-3_tsv_overlap_coef.tsv"), sep="\t")
colnames(hhsearch_OG) <- c("query_OG", "target_OG", "expect", "score", "overlap_coef")

hhsearch_OG_filter <- hhsearch_OG %>% filter(query_OG != target_OG) %>% filter(expect < 1e-3)
hhsearch_OG_filter <- hhsearch_OG_filter[order(hhsearch_OG_filter$expect, 1/hhsearch_OG_filter$score, decreasing=FALSE),]
hhsearch_OG_filter <- hhsearch_OG_filter[!duplicated(hhsearch_OG_filter[,c("query_OG", "target_OG")]),]

# Keep the top n hits per query OG
top_n_hits <- 50
hhsearch_OG_filter_ranked <- hhsearch_OG_filter %>% group_by(query_OG) %>% mutate(rank_for_query_OG = rank(expect, ties="first"))
hhsearch_OG_filter_ranked_top <- hhsearch_OG_filter_ranked[hhsearch_OG_filter_ranked$rank_for_query_OG <= top_n_hits,]
hhsearch_OG_filter <- hhsearch_OG_filter_ranked_top

# Get bidirectional hits
hhsearch_OG_filter <- hhsearch_OG_filter %>% rowwise() %>% mutate(combined_OG_ids =  paste0(sort(c(query_OG,target_OG)), collapse=","))
bidirectional_OG_ids <- hhsearch_OG_filter$combined_OG_ids[duplicated(hhsearch_OG_filter$combined_OG_ids)]
hhsearch_OG_filter_bidir <- hhsearch_OG_filter[hhsearch_OG_filter$combined_OG_ids %in% bidirectional_OG_ids,]
hhsearch_OG_filter_bidir <- hhsearch_OG_filter_bidir %>% filter(overlap_coef < 0.30)
print(paste0("Found ", nrow(hhsearch_OG_filter_bidir), " bidir hits..."))

# Get nonbidirectional hits with more stringent expect filter
hhsearch_OG_filter_nonbidir <- hhsearch_OG_filter[!hhsearch_OG_filter$combined_OG_ids %in% bidirectional_OG_ids,]
hhsearch_OG_filter_nonbidir <- hhsearch_OG_filter_nonbidir %>% filter(expect < 1e-10)
hhsearch_OG_filter_nonbidir <- hhsearch_OG_filter_nonbidir %>% filter(overlap_coef < 0.05)
print(paste0("Found ", nrow(hhsearch_OG_filter_nonbidir), " nonbidir hits..."))

# Create matrix
copies_mat <- orthogroups_copies[,2:ncol(orthogroups_copies)]
rownames(copies_mat) <- orthogroups_copies$Orthogroup

# Get subset of OGs that will potentially be merged
potential_merge_OGs <- unique(c(hhsearch_OG_filter_bidir$query_OG, hhsearch_OG_filter_bidir$target_OG, hhsearch_OG_filter_nonbidir$query_OG, hhsearch_OG_filter_nonbidir$target_OG))
copies_mat_update <- copies_mat[which(rownames(copies_mat) %in% potential_merge_OGs),]

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

# Then, merge nonbidir hits with more stringent fraction overlap and expect thresholds
print(paste0("Processing ", nrow(hhsearch_OG_filter_nonbidir), " nonbidir hits..."))
fraction_overlap_threshold <- 0.05
if (nrow(hhsearch_OG_filter_nonbidir) > 0) {
  for (i in 1:nrow(hhsearch_OG_filter_nonbidir)) {
    curr_hhsearch <- hhsearch_OG_filter_nonbidir[i,]
    
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
      
      print(paste0("Nonbidir merge (", i, " of ", nrow(hhsearch_OG_filter_nonbidir), ")"))
      
    }
  }
}

## Write out orthogroups after bidirectional and nonbidirectional merges
# saveRDS(copies_mat_update, here("data/orthogroups/first_merge", "copies_mat_update_merge.bidir.nonbidir_first_merge.rds"))

