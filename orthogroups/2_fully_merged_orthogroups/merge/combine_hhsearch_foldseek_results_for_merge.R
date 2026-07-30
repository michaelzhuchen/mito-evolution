library(tidyverse)

hhsearch_OG_filter_bidir <- read.table(here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_overlap_coef0.3_bidir_top50hits_OG_id_hhsearch_table.tsv"), sep="\t", header=TRUE)

# Read in foldseek results.
foldseek_OG <- read.table(here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_overlap_coef0.3_bidir_top50hits_mean.bitscore.prob_combined.tsv"), sep="\t")
colnames(foldseek_OG) <- c("query_OG", "target_OG", "score", "prob")

foldseek_OG$prob[which(foldseek_OG$score <= 20)] <- 0
foldseek_OG$score[is.na(foldseek_OG$score)] <- 0
foldseek_OG$score[foldseek_OG$score < 0] <- 0

hhsearch_foldseek <- merge(hhsearch_OG_filter_bidir, foldseek_OG, by.x=c("query_OG", "target_OG"), by.y=c("query_OG", "target_OG"), all.x=TRUE, all.y=TRUE)

hhsearch_foldseek$prob.y[is.na(hhsearch_foldseek$prob.y)] <- 0
hhsearch_foldseek$score.y[is.na(hhsearch_foldseek$score.y)] <- 0

# Combine probability(HMM homology | foldseek homology) using the complement
hhsearch_foldseek$prob_OR <- 1 - (1-hhsearch_foldseek$prob.x/100)*(1-hhsearch_foldseek$prob.y/100)
hhsearch_foldseek$score_sum <- hhsearch_foldseek$score.x + hhsearch_foldseek$score.y

## Write out
# write.table(hhsearch_foldseek[,c("query_OG", "target_OG", "prob_OR", "score_sum", "overlap_coef")], here("data/orthogroups/second_merge", "OG_all_merged_full.msa_vs_all_hhsuite_prob30_overlap_coef0.3_bidir_top50hits_mean.bitscore.prob_combined_processed.tsv"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


