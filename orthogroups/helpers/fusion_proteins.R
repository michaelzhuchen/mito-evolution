### Identify fusion proteins

## Read in disjoint hits for initial set of 656 euks used for raw OGs
combined_fusion_protein_filter_disjoint_original_species <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "hmmsearch_fusion_protein_candidates_expect1e-10_DBSIZE.379668_processed_1e-10_combined_sort_split_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
colnames(combined_fusion_protein_filter_disjoint_original_species) <- c("protein_id", "disjoint_hits")
combined_fusion_protein_filter_disjoint_original_species <- combined_fusion_protein_filter_disjoint_original_species[which(combined_fusion_protein_filter_disjoint_original_species$protein_id != "protein_id"),]
# Read in fusion protein candidates and OG ids
fusion_protein_candidates_OG_id <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged.hmm.foldseek_50AA_fusion_protein_candidates_OG_ids.txt"))
colnames(fusion_protein_candidates_OG_id) <- c("protein_id", "orthogroup")
# Map each protein to its OG ID
combined_fusion_protein_filter_disjoint_original_species$orthogroup <- fusion_protein_candidates_OG_id$orthogroup[match(combined_fusion_protein_filter_disjoint_original_species$protein_id, fusion_protein_candidates_OG_id$protein_id)]

## For added/updated species, singletons, mtDNA proteins
# Read in disjoint hits
combined_fusion_protein_filter_disjoint_5741 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_5741_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_81824 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_81824_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_595528 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_595528_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_667725 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_667725_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_691883 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_691883_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_946362 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_946362_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_ACANB <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_ACANB.new.proteins_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_BaSk <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_BaSk_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_CRuMs <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_CRuMs_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_parasiticplants <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_ParasiticPlants_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_mtDNA <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_mtDNA.added_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_singleton <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged.hmm.foldseek_hmm_vs_singleton_proteins_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_mtDNA_added.2025.04.07 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_mtDNA.added.2025.04.07_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_mtDNA_added.2025.05.31 <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_mtDNA.added.2025.05.31_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)
combined_fusion_protein_filter_disjoint_ltaref <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged_hmm_vs_ltaref.new.proteins_expect1e-3_DBSIZE.379668_processed_1e-10_combined.tsv_expect1e-10_disjoint.tsv"), sep="\t", header=FALSE)

# Remove incorrectly matched accessions from mtDNA and rename incorrect taxid
combined_fusion_protein_filter_disjoint_mtDNA <- combined_fusion_protein_filter_disjoint_mtDNA[!grepl("^2086695_", combined_fusion_protein_filter_disjoint_mtDNA$V1),]
combined_fusion_protein_filter_disjoint_mtDNA$V1[grep("^192875_", combined_fusion_protein_filter_disjoint_mtDNA$V1)] <- gsub("^192875_", "595528_", combined_fusion_protein_filter_disjoint_mtDNA$V1[grep("^192875_", combined_fusion_protein_filter_disjoint_mtDNA$V1)])

combined_fusion_protein_filter_disjoint_add_species <- rbind(combined_fusion_protein_filter_disjoint_5741, combined_fusion_protein_filter_disjoint_81824, combined_fusion_protein_filter_disjoint_595528, combined_fusion_protein_filter_disjoint_667725, combined_fusion_protein_filter_disjoint_691883, combined_fusion_protein_filter_disjoint_946362, combined_fusion_protein_filter_disjoint_ACANB, combined_fusion_protein_filter_disjoint_mtDNA, combined_fusion_protein_filter_disjoint_BaSk, combined_fusion_protein_filter_disjoint_CRuMs, combined_fusion_protein_filter_disjoint_parasiticplants, combined_fusion_protein_filter_disjoint_mtDNA_added.2025.04.07, combined_fusion_protein_filter_disjoint_mtDNA_added.2025.05.31, combined_fusion_protein_filter_disjoint_ltaref)
colnames(combined_fusion_protein_filter_disjoint_add_species) <- c("protein_id", "disjoint_hits")
colnames(combined_fusion_protein_filter_disjoint_singleton) <- c("protein_id", "disjoint_hits")

# Read in added fusion protein candidates with their OG ids
fusion_protein_candidates_OG_id <- read.table(here("data/orthogroups/hmmsearch_fusion_proteins", "OG_all_merged.hmm.foldseek_add.species.singleton.mtDNA_50AA_fusion_protein_candidates_OG_ids.txt"))
colnames(fusion_protein_candidates_OG_id) <- c("protein_id", "orthogroup")
# Map original member OG
combined_fusion_protein_filter_disjoint_add_species$orthogroup <- fusion_protein_candidates_OG_id$orthogroup[match(combined_fusion_protein_filter_disjoint_add_species$protein_id, fusion_protein_candidates_OG_id$protein_id)]
combined_fusion_protein_filter_disjoint_singleton$orthogroup <- fusion_protein_candidates_OG_id$orthogroup[match(combined_fusion_protein_filter_disjoint_singleton$protein_id, fusion_protein_candidates_OG_id$protein_id)]

combined_fusion_protein_filter_disjoint <- rbind(combined_fusion_protein_filter_disjoint_original_species, combined_fusion_protein_filter_disjoint_add_species, combined_fusion_protein_filter_disjoint_singleton)

combined_fusion_protein_filter_disjoint_separaterows <- combined_fusion_protein_filter_disjoint %>% separate_rows(disjoint_hits, sep=", ")

# Remove duplicates
combined_fusion_protein_filter_disjoint_separaterows <- distinct(combined_fusion_protein_filter_disjoint_separaterows)

# Filter out fusion candidates that don't have top hit to primary OG. Generally indicates a non-specific alignment, due to common domain or huge multigene family (kinase, zinc finger, very long proteins) that don't have reliable pairwise local alignments
self_hit_fusion_candidates <- combined_fusion_protein_filter_disjoint_separaterows$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows$disjoint_hits == combined_fusion_protein_filter_disjoint_separaterows$orthogroup)]
combined_fusion_protein_filter_disjoint_separaterows <- combined_fusion_protein_filter_disjoint_separaterows[combined_fusion_protein_filter_disjoint_separaterows$protein_id %in% self_hit_fusion_candidates,]

# Find new hits that aren't hit to primary OG
combined_fusion_protein_filter_disjoint_separaterows_newhit <- combined_fusion_protein_filter_disjoint_separaterows[combined_fusion_protein_filter_disjoint_separaterows$disjoint_hits != combined_fusion_protein_filter_disjoint_separaterows$orthogroup,]


## Map to updated ACANB and Tbr annotations
# Keep all ACABI transcript IDs that are identical in ACANB and are the best isoform
ACABI_fusion_newhit_protein_ids <- combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[grep("^1257118_ACABI", combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id)]
remove_ACABI_fusion_newhit_protein_ids <- ACABI_fusion_newhit_protein_ids[which(!ACABI_fusion_newhit_protein_ids %in% ACABI_ids_identical_to_ACANB_best_transcript_isoforms_accessions)]
combined_fusion_protein_filter_disjoint_separaterows_newhit <- combined_fusion_protein_filter_disjoint_separaterows_newhit[which(!combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% remove_ACABI_fusion_newhit_protein_ids),]
# Replace ACABI with ACANB protein
combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% ACANB_ACABI_transcripts$ACABI_transcript_id)] <- ACANB_ACABI_transcripts$ACANB_transcript_id[match(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% ACANB_ACABI_transcripts$ACABI_transcript_id)], ACANB_ACABI_transcripts$ACABI_transcript_id)]

# Keep all Tbr transcript ids that are best isoform in new Tbr annotation
Tbr_fusion_oldhit_protein_ids <- combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[grep("^185431_", combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id)]
remove_Tbr_fusion_oldhit_protein_ids <- Tbr_fusion_oldhit_protein_ids[which(!Tbr_fusion_oldhit_protein_ids %in% tbr_mapping$transcriptID)]
combined_fusion_protein_filter_disjoint_separaterows_newhit <- combined_fusion_protein_filter_disjoint_separaterows_newhit[which(!combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% remove_Tbr_fusion_oldhit_protein_ids),]
# Replace Tbr transcript id with gene id
combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% tbr_mapping$transcriptID)] <- tbr_mapping$ID[match(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% tbr_mapping$transcriptID)], tbr_mapping$transcriptID)]

# Keep all lta2019 that are identical in ltaref
lta_fusion_oldhit_protein_ids <- combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[grep("^5689_GET", combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id)]
remove_lta_fusion_oldhit_protein_ids <- lta_fusion_oldhit_protein_ids[which(!lta_fusion_oldhit_protein_ids %in% lta_mapping$lta2019_id)]
combined_fusion_protein_filter_disjoint_separaterows_newhit <- combined_fusion_protein_filter_disjoint_separaterows_newhit[which(!combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% remove_lta_fusion_oldhit_protein_ids),]
# Replace lta2019 id with ltaref id
combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% lta_mapping$lta2019_id)] <- lta_mapping$ltaref_id[match(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id[which(combined_fusion_protein_filter_disjoint_separaterows_newhit$protein_id %in% lta_mapping$lta2019_id)], lta_mapping$lta2019_id)]

# Get summary of new fusion genes by OG
fusion_protein_summary <- combined_fusion_protein_filter_disjoint_separaterows_newhit %>% group_by(disjoint_hits) %>% summarize(accessions = paste(unique(protein_id), sep="", collapse=", "))

