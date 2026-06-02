### Prepare training data for DeepLoc

# Load libraries
library(tidyverse)
library(Biostrings)


## Select training dataset
# Include all species
suffix <- "PCP_All_2026.04.05"
taxids <- c("5741", "185431", "5689", "3702", "1257118", "32595")

# Leave out Acanthamoeba
# suffix <- "LOO.Ac_2026.04.05"
# taxids <- c("5741", "185431", "5689", "3702", "32595")

# Leave out kinetoplastids (Trypanosoma and Leishmania)
# suffix <- "LOO.Lt.Tb_2026.04.05"
# taxids <- c("5741", "3702", "1257118", "32595")

# Leave out Arabidopsis
# suffix <- "LOO.At_2026.04.05"
# taxids <- c("5741", "185431", "5689", "1257118", "32595")

# Leave out Babesia
# suffix <- "LOO.Bd_2026.04.05"
# taxids <- c("5741", "185431", "5689", "3702", "1257118")

# Leave out Giardia
# suffix <- "LOO.Gl_2026.04.05"
# taxids <- c("185431", "5689", "3702", "1257118", "32595")


## Get protein correlation profiling (PCP) data
aca_pcp <- read.table(here('data/deeploc/mitotol', 'Ac.MitoCarta_transcript.id.txt'), header=TRUE)
ath_pcp <- read.table(here('data/deeploc/mitotol', 'At.MitoCarta.txt'), header=TRUE)
gla_pcp <- read.table(here('data/deeploc/mitotol', 'Gl.MitoCarta.txt'), header=TRUE)
lta_pcp <- read.table(here('data/deeploc/mitotol', 'Lt.MitoCarta.txt'), header=TRUE)
tbr_pcp <- read.table(here('data/deeploc/mitotol', 'Tb.MitoCarta.txt'), header=TRUE)
bdi_pcp <- read.table(here('data/deeploc/mitotol', 'Bd.MitoCarta.txt'), header=TRUE)

aca_pcp$taxid <- "1257118"
ath_pcp$taxid <- "3702"
bdi_pcp$taxid <- "32595"
gla_pcp$taxid <- "5741"
lta_pcp$taxid <- "5689"
tbr_pcp$taxid <- "185431"

aca_pcp_mito <- aca_pcp %>% filter(PCP > 0.6 | MicroscopyLabel %in% c("mito", "mito+other"))
ath_pcp_mito <- ath_pcp %>% filter(PCP > 0.6 | MicroscopyLabel %in% c("mito", "mito+other"))
bdi_pcp_mito <- bdi_pcp %>% filter(PCP > 0.96 | MicroscopyLabel %in% c("mito", "mito+other"))
gla_pcp_mito <- gla_pcp %>% filter(PCP > 0.95 | MicroscopyLabel %in% c("mito", "mito+other"))
lta_pcp_mito <- lta_pcp %>% filter(PCP > 0.6 | MicroscopyLabel %in% c("mito", "mito+other"))
tbr_pcp_mito <- tbr_pcp %>% filter(PCP > 0.6 | MicroscopyLabel %in% c("mito", "mito+other"))

aca_pcp_mito_only <- aca_pcp %>% filter(PCP > 0.95 | MicroscopyLabel %in% c("mito"))
ath_pcp_mito_only <- ath_pcp %>% filter(PCP > 0.95 | MicroscopyLabel %in% c("mito") | MitoCartaEvidence %in% c("MOM_signature"))
bdi_pcp_mito_only <- bdi_pcp %>% filter(PCP > 0.96 | MicroscopyLabel %in% c("mito"))
gla_pcp_mito_only <- gla_pcp %>% filter(PCP > 0.95 | MicroscopyLabel %in% c("mito"))
lta_pcp_mito_only <- lta_pcp %>% filter(PCP > 0.95 | MicroscopyLabel %in% c("mito"))
tbr_pcp_mito_only <- tbr_pcp %>% filter(PCP > 0.95 | MicroscopyLabel %in% c("mito"))

aca_pcp_nonmito <- aca_pcp %>% filter(PCP < 0 & PCP > -1)
ath_pcp_nonmito <- ath_pcp %>% filter(PCP < 0 & PCP > -1)
bdi_pcp_nonmito <- bdi_pcp %>% filter(PCP < 0 & PCP > -1)
gla_pcp_nonmito <- gla_pcp %>% filter(PCP < 0 & PCP > -1)
lta_pcp_nonmito <- lta_pcp %>% filter(PCP < 0 & PCP > -1)
tbr_pcp_nonmito <- tbr_pcp %>% filter(PCP < 0 & PCP > -1)

mito_pcp <- rbind(aca_pcp_mito, ath_pcp_mito, bdi_pcp_mito, gla_pcp_mito, lta_pcp_mito, tbr_pcp_mito)
mito_only_pcp <- rbind(aca_pcp_mito_only, ath_pcp_mito_only, bdi_pcp_mito_only, gla_pcp_mito_only, lta_pcp_mito_only, tbr_pcp_mito_only)
nonmito_pcp <- rbind(aca_pcp_nonmito, ath_pcp_nonmito, bdi_pcp_nonmito, gla_pcp_nonmito, lta_pcp_nonmito, tbr_pcp_nonmito)

mito_pcp <- mito_pcp %>% rowwise() %>% mutate(protein_accession = paste0(taxid, "_", ID))
mito_all_pcp <- mito_pcp
mito_only_pcp <- mito_only_pcp %>% rowwise() %>% mutate(protein_accession = paste0(taxid, "_", ID))
mito_pcp <- mito_pcp %>% filter(!protein_accession %in% mito_only_pcp$protein_accession)
nonmito_pcp <- nonmito_pcp %>% rowwise() %>% mutate(protein_accession = paste0(taxid, "_", ID))

# Exclude any MitoCarta proteins from nonmito set.
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2026.04.05.tsv"), sep="\t", header=TRUE)
nonmito_pcp <- nonmito_pcp %>% filter(!protein_accession %in% gold_gene_accession_OG_id_df$gene_accession)


## Process ggsearch results to add mitoepi proteins to original swissprot training dataset

# Read in Swissprot labels from original training data
multisub_idmapping <- read.delim(here("data/deeploc/swissprot", "multisub_5_partitions_unique_idmapping_2025_04_30.tsv"), header=TRUE)
multisub_labels <- read.csv(here("data/deeploc/swissprot", "multisub_5_partitions_unique.csv"))
multisub_labels$species <- multisub_idmapping$Organism[match(multisub_labels$ACC, multisub_idmapping$From)]
multisub_labels$taxid <- multisub_idmapping$Organism..ID.[match(multisub_labels$ACC, multisub_idmapping$From)]

# Read in sequences
swissprot_fasta <- readAAStringSet(here("data/deeploc/swissprot", "deeploc_swissprot.fasta"))
mitoepi_fasta <- readAAStringSet(here("data/deeploc/mitotol", "mitoepi_species_combined.fasta")) # all MITO-EPI species proteins
combined_fasta <- c(mitoepi_fasta, swissprot_fasta)
names(combined_fasta)[grep("^1257118_", names(combined_fasta))] <- gsub("_t.*", "", names(combined_fasta)[grep("^1257118_", names(combined_fasta))])

## Update human + yeast mito/nonmito labels in preexisting Swissprot data
# Ensure that human mito/nonmito labels are correct
multisub_labels_human <- multisub_labels %>% filter(taxid == "9606")
human_mito_uniprot <- read.delim(here("data/orthogroups/idmapping", "human_mitocarta_geneid_to_uniprot_idmapping_2024_10_01.tsv"))
multisub_labels_human$Mitochondrion[multisub_labels_human$ACC %in% human_mito_uniprot$Entry] <- 1
multisub_labels_human$Mitochondrion[!multisub_labels_human$ACC %in% human_mito_uniprot$Entry] <- 0

# Ensure that yeast mito/nonmito labels are correct
multisub_labels_yeast <- multisub_labels %>% filter(taxid == "559292")
yeast_mito_uniprot <- read.delim(here("data/orthogroups/idmapping", "sce_mitocarta_genename_to_uniprot_idmapping_2024_10_01.tsv"))
multisub_labels_yeast$Mitochondrion[multisub_labels_yeast$ACC %in% yeast_mito_uniprot$Entry] <- 1
multisub_labels_yeast$Mitochondrion[!multisub_labels_yeast$ACC %in% yeast_mito_uniprot$Entry] <- 0

multisub_labels_nohumanyeast <- multisub_labels %>% filter(!taxid %in% c("9606", "559292"))
multisub_labels <- rbind(multisub_labels_nohumanyeast, multisub_labels_human, multisub_labels_yeast)

## Add missing human and yeast sequences to preexisting Swissprot data
# Add missing human mito
human_fasta <- readAAStringSet(here("species_fastas", "9606.fasta"))
human_mito_uniprot <- human_mito_uniprot %>% mutate(accession = paste0("9606_", Entry.Name))
human_mito_uniprot_new <- human_mito_uniprot %>% filter(!Entry %in% multisub_labels_human$ACC & !accession %in% multisub_labels_human$ACC) %>% filter(accession %in% names(human_fasta))
new_protein_accessions <- unique(human_mito_uniprot_new$accession)
multisub_labels_mito_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_human))
colnames(multisub_labels_mito_new) <- colnames(multisub_labels_human)
multisub_labels_mito_new <- as.data.frame(multisub_labels_mito_new)
multisub_labels_mito_new$ACC <- new_protein_accessions
if (nrow(multisub_labels_mito_new) > 0) {
  multisub_labels_mito_new$Mitochondrion <- 1
  multisub_labels_mito_new$taxid <- "9606"
  multisub_labels_mito_new$Sequence <- as.character(combined_fasta[multisub_labels_mito_new$ACC])
}
multisub_labels <- rbind(multisub_labels, multisub_labels_mito_new)
# Add missing human nonmito
human_uniprot_all <- read.delim(here("data/orthogroups/idmapping", "human_complete_geneid_to_uniprot_idmapping_2024_10_25.tsv"), header=TRUE)
human_uniprot_all <- human_uniprot_all %>% mutate(accession = paste0("9606_", Entry.Name))
multisub_labels_human <- multisub_labels %>% filter(taxid == "9606")
human_uniprot_all_new <- human_uniprot_all %>% filter(!Entry %in% multisub_labels_human$ACC & !accession %in% multisub_labels_human$ACC) %>% filter(accession %in% names(human_fasta))
new_protein_accessions <- unique(human_uniprot_all_new$accession)
multisub_labels_nonmito_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_human))
colnames(multisub_labels_nonmito_new) <- colnames(multisub_labels_human)
multisub_labels_nonmito_new <- as.data.frame(multisub_labels_nonmito_new)
multisub_labels_nonmito_new$ACC <- new_protein_accessions
if (nrow(multisub_labels_nonmito_new) > 0) {
  multisub_labels_nonmito_new$Mitochondrion <- 0
  multisub_labels_nonmito_new$taxid <- "9606"
  multisub_labels_nonmito_new$Sequence <- as.character(combined_fasta[multisub_labels_nonmito_new$ACC])
}
multisub_labels <- rbind(multisub_labels, multisub_labels_nonmito_new)

# Add missing yeast mito
yeast_fasta <- readAAStringSet(here("species_fastas", "559292.fasta"))
yeast_mito_uniprot <- yeast_mito_uniprot %>% mutate(accession = paste0("559292_", Entry.Name))
yeast_mito_uniprot_new <- yeast_mito_uniprot %>% filter(!Entry %in% multisub_labels_yeast$ACC & !accession %in% multisub_labels_yeast$ACC) %>% filter(accession %in% names(yeast_fasta))
new_protein_accessions <- unique(yeast_mito_uniprot_new$accession)
multisub_labels_mito_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_yeast))
colnames(multisub_labels_mito_new) <- colnames(multisub_labels_yeast)
multisub_labels_mito_new <- as.data.frame(multisub_labels_mito_new)
multisub_labels_mito_new$ACC <- new_protein_accessions
if (nrow(multisub_labels_mito_new) > 0) {
  multisub_labels_mito_new$Mitochondrion <- 1
  multisub_labels_mito_new$taxid <- "559292"
  multisub_labels_mito_new$Sequence <- as.character(combined_fasta[multisub_labels_mito_new$ACC])
}
multisub_labels <- rbind(multisub_labels, multisub_labels_mito_new)
# Add missing yeast nonmito
yeast_uniprot_all <- read.delim(here("data/orthogroups/idmapping", "sce_geneid_to_uniprot_idmapping_2024_10_25.tsv"), header=TRUE)
yeast_uniprot_all$accession <- paste0("559292_", yeast_uniprot_all$Entry.Name)
yeast_uniprot_all <- yeast_uniprot_all %>% filter(accession %in% names(yeast_fasta))
multisub_labels_yeast <- multisub_labels %>% filter(taxid == "559292")
yeast_uniprot_all_new <- yeast_uniprot_all %>% filter(!Entry %in% multisub_labels_yeast$ACC & !accession %in% multisub_labels_yeast$ACC)
new_protein_accessions <- unique(yeast_uniprot_all_new$accession)
multisub_labels_nonmito_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_yeast))
colnames(multisub_labels_nonmito_new) <- colnames(multisub_labels_yeast)
multisub_labels_nonmito_new <- as.data.frame(multisub_labels_nonmito_new)
multisub_labels_nonmito_new$ACC <- new_protein_accessions
if (nrow(multisub_labels_nonmito_new) > 0) {
  multisub_labels_nonmito_new$Mitochondrion <- 0
  multisub_labels_nonmito_new$taxid <- "559292"
  multisub_labels_nonmito_new$Sequence <- as.character(combined_fasta[multisub_labels_nonmito_new$ACC])
}
multisub_labels <- rbind(multisub_labels, multisub_labels_nonmito_new)

# Add mito-epi species proteins
pident_threshold <- 90
pident_partition_threshold <- 30

new_multisub_labels <- c()
old_accessions_to_replace <- c()

ggsearch <- read.table(here("data/deeploc/ggsearch", "ggsearch_swissprot_mitoepi_combined_fasta_vs_self_expect1_combined.out"), sep="\t", header=FALSE)
colnames(ggsearch) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

for (curr_taxid in taxids) {
  # Get identical/near identical top hits for current species
  ggsearch_filter <- ggsearch %>% filter(gsub("_.*", "", qseqid) == curr_taxid) %>% filter(evalue < 1e-3) %>% filter(pident >= pident_threshold) %>% group_by(qseqid) %>% filter(evalue == min(evalue) & bitscore == max(bitscore))
  ggsearch_filter$target_species <- multisub_idmapping$Organism[match(ggsearch_filter$sseqid, multisub_idmapping$From)]
  
  # Mito
  ggsearch_filter_mito <- ggsearch_filter %>% filter(qseqid %in% mito_pcp$protein_accession)
  multisub_labels_mito <- multisub_labels %>% filter(ACC %in% ggsearch_filter_mito$sseqid) %>% filter(taxid == curr_taxid)
  old_accessions_to_replace <- c(old_accessions_to_replace, multisub_labels_mito$ACC)
  multisub_labels_mito$ACC <- ggsearch_filter_mito$qseqid[match(multisub_labels_mito$ACC, ggsearch_filter_mito$sseqid)]
  
  mito_pcp_curr <- mito_pcp %>% filter(taxid == curr_taxid)
  mito_pcp_curr_new <- mito_pcp_curr %>% filter(!protein_accession %in% multisub_labels_mito$ACC)
  new_protein_accessions <- unique(mito_pcp_curr_new$protein_accession)
  multisub_labels_mito_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_mito))
  colnames(multisub_labels_mito_new) <- colnames(multisub_labels_mito)
  multisub_labels_mito_new <- as.data.frame(multisub_labels_mito_new)
  multisub_labels_mito_new$ACC <- new_protein_accessions
  multisub_labels_mito <- rbind(multisub_labels_mito, multisub_labels_mito_new)
  
  if (nrow(multisub_labels_mito) > 0) {
    multisub_labels_mito$Mitochondrion <- 1
  }
  
  # Nonmito
  ggsearch_filter_nonmito <- ggsearch_filter %>% filter(qseqid %in% nonmito_pcp$protein_accession)
  multisub_labels_nonmito <- multisub_labels %>% filter(ACC %in% ggsearch_filter_nonmito$sseqid) %>% filter(taxid == curr_taxid)
  old_accessions_to_replace <- c(old_accessions_to_replace, multisub_labels_nonmito$ACC)
  multisub_labels_nonmito$ACC <- ggsearch_filter_nonmito$qseqid[match(multisub_labels_nonmito$ACC, ggsearch_filter_nonmito$sseqid)]
  
  nonmito_pcp_curr <- nonmito_pcp %>% filter(taxid == curr_taxid)
  nonmito_pcp_curr_new <- nonmito_pcp_curr %>% filter(!protein_accession %in% multisub_labels_nonmito$ACC)
  new_protein_accessions <- unique(nonmito_pcp_curr_new$protein_accession)
  multisub_labels_nonmito_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_nonmito))
  colnames(multisub_labels_nonmito_new) <- colnames(multisub_labels_nonmito)
  multisub_labels_nonmito_new <- as.data.frame(multisub_labels_nonmito_new)
  multisub_labels_nonmito_new$ACC <- new_protein_accessions
  multisub_labels_nonmito <- rbind(multisub_labels_nonmito, multisub_labels_nonmito_new)
  
  if (nrow(multisub_labels_nonmito) > 0) {
    multisub_labels_nonmito$Mitochondrion <- 0
  }
  
  # Mito Only
  ggsearch_filter_mito_only <- ggsearch_filter %>% filter(qseqid %in% mito_only_pcp$protein_accession)
  multisub_labels_mito_only <- multisub_labels %>% filter(ACC %in% ggsearch_filter_mito_only$sseqid) %>% filter(taxid == curr_taxid)
  old_accessions_to_replace <- c(old_accessions_to_replace, multisub_labels_mito_only$ACC)
  multisub_labels_mito_only$ACC <- ggsearch_filter_mito_only$qseqid[match(multisub_labels_mito_only$ACC, ggsearch_filter_mito_only$sseqid)]
  
  mito_only_pcp_curr <- mito_only_pcp %>% filter(taxid == curr_taxid)
  mito_only_pcp_curr_new <- mito_only_pcp_curr %>% filter(!protein_accession %in% multisub_labels_mito_only$ACC)
  new_protein_accessions <- unique(mito_only_pcp_curr_new$protein_accession)
  multisub_labels_mito_only_new <- matrix(data = NA, nrow = length(new_protein_accessions), ncol = ncol(multisub_labels_mito_only))
  colnames(multisub_labels_mito_only_new) <- colnames(multisub_labels_mito_only)
  multisub_labels_mito_only_new <- as.data.frame(multisub_labels_mito_only_new)
  multisub_labels_mito_only_new$ACC <- new_protein_accessions
  multisub_labels_mito_only <- rbind(multisub_labels_mito_only, multisub_labels_mito_only_new)
  
  multisub_labels_mito_only[,6:15] <- 0 # Assign all non-mito labels to 0
  multisub_labels_mito_only$Mitochondrion <- 1
  
  
  new_multisub_labels_curr <- rbind(multisub_labels_mito, multisub_labels_mito_only, multisub_labels_nonmito)
  
  new_multisub_labels_curr$taxid <- curr_taxid
  new_multisub_labels_curr$species <- NA
  
  # Fetch sequences
  new_multisub_labels_curr <- new_multisub_labels_curr %>% filter(ACC %in% names(combined_fasta))
  new_multisub_labels_curr$Sequence <- as.character(combined_fasta[new_multisub_labels_curr$ACC])
  
  
  new_multisub_labels <- rbind(new_multisub_labels, new_multisub_labels_curr)
}

multisub_labels_swissprot_retained <- multisub_labels %>% filter(!ACC %in% old_accessions_to_replace)
multisub_labels_swissprot_mitoepi <- rbind(multisub_labels_swissprot_retained, new_multisub_labels)

# Exclude organelle DNA encoded proteins
mtdna_proteins <- read.table(here("data/deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1
nonmito_organelle_proteins <- read.table(here("data/deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
cpdna_proteins <- grep("ATCG", new_multisub_labels$ACC, value=TRUE)
organelle_protein_accessions <- unique(c(mtdna_proteins, nonmito_organelle_proteins, cpdna_proteins))
new_multisub_labels <- new_multisub_labels %>% filter(!ACC %in% organelle_protein_accessions)
multisub_labels_swissprot_retained <- multisub_labels_swissprot_retained %>% filter(!ACC %in% organelle_protein_accessions)
multisub_labels_swissprot_mitoepi <- multisub_labels_swissprot_mitoepi %>% filter(!ACC %in% organelle_protein_accessions)

# Keep sequences >= 40 AA
new_multisub_labels <- new_multisub_labels %>% rowwise() %>% filter(nchar(Sequence) >= 40)
multisub_labels_swissprot_retained <- multisub_labels_swissprot_retained %>% rowwise() %>% filter(nchar(Sequence) >= 40)
multisub_labels_swissprot_mitoepi <- multisub_labels_swissprot_mitoepi %>% rowwise() %>% filter(nchar(Sequence) >= 40)

# Remove duplicates
new_multisub_labels <- new_multisub_labels %>% distinct()
multisub_labels_swissprot_retained <- multisub_labels_swissprot_retained %>% distinct()
multisub_labels_swissprot_mitoepi <- multisub_labels_swissprot_mitoepi %>% distinct()

# Sanity check
multisub_labels_swissprot_mitoepi %>% filter(taxid %in% c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595")) %>% group_by(taxid) %>% summarize(n_mito = sum(Mitochondrion == 1), n_nonmito = sum(Mitochondrion == 0), n_total = n())


# Prepare header-labeled fasta for graph-part
combined_fasta_for_train <- combined_fasta[multisub_labels_swissprot_mitoepi$ACC]

# Assign labels, in order of most frequent to least frequent, so that least frequent labels are prioritized by graph-part (which only considers the last label in the header). Add Mitochondrion label last to prioritize.
sort(colSums(multisub_labels_swissprot_mitoepi[,6:15], na.rm=TRUE), decreasing=TRUE)
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Cytoplasm == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Cytoplasm == 1)], "|label=Cytoplasm")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Nucleus == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Nucleus == 1)], "|label=Nucleus")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Cell.membrane == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Cell.membrane == 1)], "|label=Cell.membrane")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Extracellular == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Extracellular == 1)], "|label=Extracellular")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Endoplasmic.reticulum == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Endoplasmic.reticulum == 1)], "|label=Endoplasmic.reticulum")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Lysosome.Vacuole == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Lysosome.Vacuole == 1)], "|label=Lysosome.Vacuole")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Golgi.apparatus == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Golgi.apparatus == 1)], "|label=Golgi.apparatus")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Plastid == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Plastid == 1)], "|label=Plastid")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Peroxisome == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Peroxisome == 1)], "|label=Peroxisome")
names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Mitochondrion == 1)] <- paste0(names(combined_fasta_for_train)[which(multisub_labels_swissprot_mitoepi$Mitochondrion == 1)], "|label=Mitochondrion")


## Write out
# write.csv(multisub_labels_swissprot_mitoepi, here("data/deeploc/data_files", paste0("multisub_swissprot.retained_added.", suffix, ".csv")), row.names=FALSE, quote=FALSE)
# writeXStringSet(combined_fasta_for_train, here("data/deeploc/data_files", paste0("swissprot_", suffix, ".fasta")))

