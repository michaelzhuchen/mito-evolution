### Get added proteins from added/updated species and added mtDNA annotations

# Added/updated species
tophit_5741 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_5741_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_81824 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_81824_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_595528 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_595528_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_667725 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_667725_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_691883 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_691883_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_946362 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_946362_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_ACANB_new_proteins <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_ACANB.new.proteins_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_BaSk <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_BaSk_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_CRuMs <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_CRuMs_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_parasiticplants <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_ParasiticPlants_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_ltaref_new_proteins <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_ltaref.new.proteins_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))

# Added mtDNA annotations
tophit_mtDNA_new_proteins <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_mtDNA.added_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_mtDNA_new_proteins_2025.04.07 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_mtDNA.added.2025.04.07_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))
tophit_mtDNA_new_proteins_2025.05.31 <- read.table(here("data/orthogroups/hmmsearch_added_proteins", "OG_all_merged_hmm_vs_mtDNA.added.2025.05.31_expect1e-3_DBSIZE.379668_combined.out_tophit.tsv"))

# Remove incorrectly matched accessions from mtDNA and rename incorrect taxid
tophit_mtDNA_new_proteins <- tophit_mtDNA_new_proteins[!grepl("^2086695_", tophit_mtDNA_new_proteins$V1),]
tophit_mtDNA_new_proteins$V1[grep("^192875_", tophit_mtDNA_new_proteins$V1)] <- gsub("^192875_", "595528_", tophit_mtDNA_new_proteins$V1[grep("^192875_", tophit_mtDNA_new_proteins$V1)])

tophit_add_species <- rbind(tophit_5741, tophit_81824, tophit_595528, tophit_667725, tophit_691883, tophit_946362, tophit_ACANB_new_proteins, tophit_mtDNA_new_proteins, tophit_BaSk, tophit_CRuMs, tophit_parasiticplants, tophit_mtDNA_new_proteins_2025.04.07, tophit_ltaref_new_proteins, tophit_mtDNA_new_proteins_2025.05.31)
colnames(tophit_add_species) <- c("protein_id", "OG_id")

# Split ties into separate rows
tophit_add_species_separaterows <- tophit_add_species %>% separate_rows(OG_id, sep=",")

# Get summary of new proteins by OG to add to their new OGs in block below
tophit_add_species_summary <- tophit_add_species_separaterows %>% group_by(OG_id) %>% summarize(accessions = paste(unique(protein_id), sep="", collapse=", "))



