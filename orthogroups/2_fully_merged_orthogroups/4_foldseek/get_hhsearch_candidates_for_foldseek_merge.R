library(tidyverse)

# Read in OG_all_merged hhsearch results with precomputed overlap coef using the updated OG all merged copies mat.
hhsearch_OG <- read.table(here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_tsv_overlap_coef_with_merged_combined.tsv"), sep="\t")
colnames(hhsearch_OG) <- c("query_OG", "target_OG", "expect", "score", "overlap_coef", "prob")

# Filter
hhsearch_OG_filter <- hhsearch_OG %>% filter(query_OG != target_OG) %>% filter(overlap_coef < 0.3)
hhsearch_OG_filter <- hhsearch_OG_filter[order(hhsearch_OG_filter$expect, 1/hhsearch_OG_filter$score, decreasing=FALSE),]
hhsearch_OG_filter <- hhsearch_OG_filter[!duplicated(hhsearch_OG_filter[,c("query_OG", "target_OG")]),]

# Keep the top n hits per query OG
top_n_hits <- 50
hhsearch_OG_filter_ranked <- hhsearch_OG_filter %>% group_by(query_OG) %>% mutate(rank_for_query_OG = rank(expect, ties="first"))
hhsearch_OG_filter_ranked_top <- hhsearch_OG_filter_ranked[hhsearch_OG_filter_ranked$rank_for_query_OG <= top_n_hits,]
hhsearch_OG_filter <- hhsearch_OG_filter_ranked_top

# Identify unique pairs
hhsearch_OG_filter <- hhsearch_OG_filter %>% rowwise() %>% mutate(combined_OG_ids = paste0(sort(c(query_OG,target_OG)), collapse=","))

# # Get bidirectional hits
bidirectional_OG_ids <- hhsearch_OG_filter$combined_OG_ids[duplicated(hhsearch_OG_filter$combined_OG_ids)]

hhsearch_OG_filter_bidir <- hhsearch_OG_filter[hhsearch_OG_filter$combined_OG_ids %in% bidirectional_OG_ids,]

# Remove duplicates, taking the top ranked hit
hhsearch_OG_filter_bidir <- hhsearch_OG_filter_bidir[!duplicated(hhsearch_OG_filter_bidir$combined_OG_ids),]

## Write out filtered OG pairs for foldseek
# write.table(hhsearch_OG_filter_bidir[, c("query_OG", "target_OG")], here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_overlap_coef0.3_bidir_top50hits_OG_id_query.target_v2.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# write.table(hhsearch_OG_filter_bidir, here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_overlap_coef0.3_bidir_top50hits_OG_id_hhsearch_table.tsv"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

