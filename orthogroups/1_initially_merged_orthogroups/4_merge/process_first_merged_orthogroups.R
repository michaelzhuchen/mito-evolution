## Process merged orthogroups

# Read in raw OGs
orthogroups <- read.delim(here("data/orthogroups/raw_orthogroups/Orthogroups", "Orthogroups.tsv"), header=TRUE)
all_OG_ids <- read.table(here("data/orthogroups/raw_orthogroups/Orthogroups", "raw_OG_ids_652euks.txt"))$V1
# Exclude the OGs that are only present in the excluded species
orthogroups <- orthogroups[which(orthogroups$Orthogroup %in% all_OG_ids),]

# Read in merged copies mat
copies_mat_update <- readRDS(here("data/orthogroups/first_merge", "copies_mat_update_merge.bidir.nonbidir_first_merge.rds"))

# Remove the original OG rows that have been merged. These start with an "X"
copies_mat_update <- copies_mat_update[!grepl("X", rownames(copies_mat_update)),]

merged_OG_ids <- rownames(copies_mat_update)[grep(",", rownames(copies_mat_update))]
merged_OG_df <- data.frame(OG_index=1:length(merged_OG_ids), OG_id=merged_OG_ids)
merged_OG_seprows <- merged_OG_df %>% separate_rows(OG_id, sep=",")

# Retrieve the raw OGs to be merged
orthogroups_merge <- orthogroups[which(orthogroups$Orthogroup %in% merged_OG_seprows$OG_id),]
orthogroups_merge$merge_index <- merged_OG_seprows$OG_index[match(orthogroups_merge$Orthogroup, merged_OG_seprows$OG_id)]
orthogroups_merge <- orthogroups_merge %>% group_by(merge_index) %>% summarize(across(everything(), ~paste(.[. != ""], collapse = ", ")))
orthogroups_merge <- orthogroups_merge[,2:ncol(orthogroups_merge)]

## Write out
# write.table(orthogroups_merge, here("data/orthogroups/first_merge", "first_merged_orthogroups.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


