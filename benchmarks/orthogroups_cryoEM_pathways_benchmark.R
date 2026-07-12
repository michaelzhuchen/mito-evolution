library(here)
library(tidyverse)
library(ggplot2)
library(reshape2)

# Set ggplot theme
theme_set(theme_classic())

## Specify whether to use PhROGs (+/- abSENSE-HMM predictions)
BOOL_USE_PHROG <- TRUE
BOOL_USE_ABSENSE <- FALSE

# Read in mappings from Uniprot ID to Swissprot entry
swissprot_mapping <- read.table(here("data", "orthogroups", "idmapping", "uniprot_swissprot_id_mapping.txt"), sep=" ")
colnames(swissprot_mapping) <- c("uniprot_entry", "uniprot_accession")
# Remove duplicates to keep just the current swissprot main accession
swissprot_mapping <- swissprot_mapping[!duplicated(swissprot_mapping$uniprot_entry),]

## Read in datasets
ath_idmapping <- read.delim(here("data", "orthogroups", "idmapping", "ath_idmapping_2024_11_08.tsv"))

# Read in new Tbr mapping for fasta with best isoform per gene
tbr_mapping <- read.table(here("data", "orthogroups", "idmapping", "tbrgene.mapping.txt"), header=TRUE)
tbr_mapping$transcriptID <- gsub(":", "_", tbr_mapping$transcriptID, fixed=TRUE)
tbr_mapping$ID <- paste0("185431_", tbr_mapping$ID)
tbr_mapping$transcriptID <- paste0("185431_", tbr_mapping$transcriptID)

# Refined OGs mapped to uniprot
orthogroups_collapse_accessions_seprows_uniprot <- read.table(here("data", "orthogroups", "idmapping", "refined_OGs_euk673spp_uniprot_mapping_wide.txt"), sep="\t", header=FALSE)
colnames(orthogroups_collapse_accessions_seprows_uniprot) <- c("OG_id", "accessions", "uniprot_accession")
# Map tbr to tbrgene
orthogroups_collapse_accessions_seprows_uniprot$accessions[orthogroups_collapse_accessions_seprows_uniprot$accessions %in% tbr_mapping$transcriptID] <- tbr_mapping$ID[match(orthogroups_collapse_accessions_seprows_uniprot$accessions[orthogroups_collapse_accessions_seprows_uniprot$accessions %in% tbr_mapping$transcriptID], tbr_mapping$transcriptID)]

# Read in blast data
blast_out_raw <- read.table(here("data", "benchmarks", "orthogroups", "blast_human_vs_keyspecies_forward.reverse.out"), sep="\t", header=FALSE)
colnames(blast_out_raw) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

blast_out_raw$query_taxid <- gsub("_.*", "", blast_out_raw$qseqid)
blast_out_raw$target_taxid <- gsub("_.*", "", blast_out_raw$sseqid)

# Keep top hits by expect, break ties with bitscore
blast_out_raw <- blast_out_raw %>% group_by(qseqid, target_taxid) %>% filter(evalue == min(evalue)) %>% filter(bitscore == max(bitscore))

key_species <- c("9606", "559292", "185431", "3702", "5741", "3562", "508771", "312017")
blast_out_raw <- blast_out_raw[blast_out_raw$query_taxid %in% key_species,]
blast_out_raw <- blast_out_raw[blast_out_raw$target_taxid %in% key_species,]

blast_out_raw$query_uniprot_entry <- gsub("9606_", "", blast_out_raw$qseqid, fixed=TRUE)
blast_out_raw$query_uniprot_accession <- blast_out_raw$query_uniprot_entry
blast_out_raw$target_uniprot_entry <- sub("^[^_]*_", "", blast_out_raw$sseqid)
blast_out_raw$target_uniprot_accession <- blast_out_raw$target_uniprot_entry

blast_out_tophit <- blast_out_raw %>% filter(evalue < 1e-3)

# Find BBHs with more lenient expect
blast_out_BBH <- blast_out_raw[which(blast_out_raw$evalue < 1),]
blast_out_BBH <- blast_out_BBH %>% rowwise() %>% mutate(query_target_pair = paste0(sort(c(qseqid, sseqid)), collapse=","))
BBH_query_target_pairs <- unique(blast_out_BBH$query_target_pair[duplicated(blast_out_BBH$query_target_pair)])
# Find BBHs by query target pairs or if it's a human sequence with identical qseqid and sseqid
blast_out_BBH <- blast_out_BBH[blast_out_BBH$query_target_pair %in% BBH_query_target_pairs | blast_out_BBH$qseqid == blast_out_BBH$sseqid,]
# Remove duplicates
blast_out_BBH <- blast_out_BBH[!duplicated(blast_out_BBH[,c("qseqid", "sseqid")]),]


# Read in PLMSearch+PLMAlign data
plmsearch_out_raw <- read.table(here("data", "benchmarks", "orthogroups", "plmalign_score_sort_forward_reverse_combined.txt"), sep="\t", header=FALSE)
colnames(plmsearch_out_raw) <- c("query", "target", "score")
plmsearch_out_raw <- plmsearch_out_raw %>% mutate(query_taxid = gsub("_.*", "", query), target_taxid = gsub("_.*", "", target))

plmsearch_out_raw$query_uniprot_entry <- gsub("9606_", "", plmsearch_out_raw$query, fixed=TRUE)
plmsearch_out_raw$query_uniprot_accession <- plmsearch_out_raw$query_uniprot_entry
plmsearch_out_raw$target_uniprot_entry <- sub("^[^_]*_", "", plmsearch_out_raw$target)
plmsearch_out_raw$target_uniprot_accession <- plmsearch_out_raw$target_uniprot_entry

plmsearch_out_tophit <- plmsearch_out_raw %>% group_by(query, target_taxid) %>% filter(score == max(score))

# Find BBHs with more lenient expect
plmsearch_out_BBH <- plmsearch_out_tophit
plmsearch_out_BBH <- plmsearch_out_BBH %>% rowwise() %>% mutate(query_target_pair = paste0(sort(c(query, target)), collapse=","))
BBH_query_target_pairs <- unique(plmsearch_out_BBH$query_target_pair[duplicated(plmsearch_out_BBH$query_target_pair)])
# Find BBHs by query target pairs or if it's a human sequence with identical qseqid and sseqid
plmsearch_out_BBH <- plmsearch_out_BBH[plmsearch_out_BBH$query_target_pair %in% BBH_query_target_pairs | plmsearch_out_BBH$query == plmsearch_out_BBH$target,]
# Remove duplicates
plmsearch_out_BBH <- plmsearch_out_BBH[!duplicated(plmsearch_out_BBH[,c("query", "target")]),]


# Read in Foldseek data
foldseek_out_raw <- read.table(here("data", "benchmarks", "orthogroups", "human_vs_model_species_proteomes_AFDB_v6_forward_reverse_combined_reduced.tsv"), sep="\t", header=TRUE)

# Keep top hits by expect, break ties with bitscore
foldseek_out_raw <- foldseek_out_raw %>% group_by(query, target_taxid) %>% filter(evalue == min(evalue)) %>% filter(bits == max(bits))

foldseek_out_tophit <- foldseek_out_raw %>% filter(evalue < 1e-3)

# Find BBHs with more lenient expect
foldseek_out_BBH <- foldseek_out_raw[which(foldseek_out_raw$evalue < 1),]
foldseek_out_BBH <- foldseek_out_BBH %>% rowwise() %>% mutate(query_target_pair = paste0(sort(c(query, target)), collapse=","))
BBH_query_target_pairs <- unique(foldseek_out_BBH$query_target_pair[duplicated(foldseek_out_BBH$query_target_pair)])
# Find BBHs by query target pairs or if it's a human sequence with identical qseqid and sseqid
foldseek_out_BBH <- foldseek_out_BBH[foldseek_out_BBH$query_target_pair %in% BBH_query_target_pairs | foldseek_out_BBH$query == foldseek_out_BBH$target,]
# Remove duplicates
foldseek_out_BBH <- foldseek_out_BBH[!duplicated(foldseek_out_BBH[,c("query", "target")]),]



# Read in EggNOG6 data at Eukaryota level
root_OG <- read.table(here("data", "benchmarks", "orthogroups", "eggnog6_euks.e6.og2seqs_and_species.tsv"), sep="\t")
colnames(root_OG) <- c("tax_id", "OG_id", "num_species", "num_proteins", "species_ids", "protein_ids")
eggnog_reference_species <- c("9606", "4932", "3702", "185431", "184922", "3562", "312017")
root_OG_long <- root_OG %>% separate_rows(protein_ids, sep=",")
root_OG_long <- root_OG_long %>% filter(gsub("\\..*", "", protein_ids) %in% eggnog_reference_species)


# Read in OrthoDB data at Eukaryota level
odb_genes <- read.delim(here("data", "benchmarks", "orthogroups", "odb12v1_genes_keyspecies.tab"), header=FALSE)
colnames(odb_genes)[1] <- c("odb_gene_id")
colnames(odb_genes)[2] <- c("species_id")
colnames(odb_genes)[3] <- c("seq_id")
colnames(odb_genes)[4] <- c("synonym")
colnames(odb_genes)[5] <- c("uniprot_id")
colnames(odb_genes)[6] <- c("gene_id")
odb_genes$gene_id <- gsub("\\..*", "", odb_genes$gene_id)
odb_eukaryota_ogs <- read.delim(here("data", "benchmarks", "orthogroups", "odb12v1_OG2genes_at2759_keyspecies.tab"), header=FALSE)
colnames(odb_eukaryota_ogs) <- c("OG_id", "odb_gene_id")
odb_eukaryota_ogs_genes <- merge(odb_eukaryota_ogs, odb_genes, by="odb_gene_id", all=TRUE)
odb_eukaryota_ogs_genes$OG_id[is.na(odb_eukaryota_ogs_genes$OG_id)] <- 1:sum(is.na(odb_eukaryota_ogs_genes$OG_id))
# OrthoDB ID mapping for Trypanosoma
idmapping_tbr <- read.delim(here("data", "orthogroups", "idmapping", "tbr_idmapping_2024_11_08.tsv"), header=TRUE)
idmapping_tbr$From <- gsub("TriTrypDB:", "", idmapping_tbr$From, fixed=TRUE)
# OrthoDB ID mapping for Arabidopsis
idmapping_ath <- read.delim(here("data", "orthogroups", "idmapping", "ath_uniprot_to_uniprot_idmapping_2025_07_02.tsv"), header=TRUE)
idmapping_ath <- idmapping_ath %>% separate_rows(TAIR, sep = ";") %>% separate_rows(Araport, sep = ";") %>% separate_rows(KEGG, sep = ";")
idmapping_ath$KEGG <- gsub("ath:", "", idmapping_ath$KEGG, fixed=TRUE)

# Read in PhROGs
phrogs_long <- read.table(here("data", "phylogenetically_resolved_orthogroups", "species_tree_1", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long <- phrogs_long %>% mutate(OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id)) %>% ungroup()

# Read in absense data
homology_power_agg <- read.table(here("data", "abSENSE_HMM", "species_tree_1", "absense_results.tsv"), sep="\t", header=TRUE)
colnames(homology_power_agg) <- c("PhROG_id", "human_accessions", "L", "R", "OG_mrca_label", "OG_outgroup_mrca_label", "outgroup_kingdom_species", "fraction_of_detectable_outgroup_kingdoms", "probability_of_detection_in_any_outgroup_species", "hmm_length", "r_squared", "comments")
homology_power_per_species <- read.table(here("data", "abSENSE_HMM", "species_tree_1", "absense_results_per_species.tsv"), sep="\t", header=FALSE)
colnames(homology_power_per_species) <- c("PhROG_id", "taxid", "probability_of_detection")
homology_power_agg <- homology_power_agg %>% rowwise() %>% mutate(OG_id = gsub("_.*", "", PhROG_id))
homology_power_agg_seprows <- homology_power_agg %>% separate_rows(human_accessions, sep=",")


## MitoRibo
# Read in mitoribo gene sets
mitoribo_table <- read.csv(here("data", "annotation", "pathways", "mitoribosome_structure_annotation.csv"), header=TRUE)

mitoribo_table_core_bacterial <- mitoribo_table[which(mitoribo_table$type == "core bacterial" & mitoribo_table$sum == 4),]
mitoribo_table_shared_euk <- mitoribo_table[which(mitoribo_table$type == "shared euk" & mitoribo_table$sum == 4),]

mitoribo_table_core <- rbind(mitoribo_table_core_bacterial, mitoribo_table_shared_euk)


# OG data mitoribo
OG_mitoribo <- c()
for (i in 1:nrow(mitoribo_table_core)) {
  mitoribo_curr_row <- mitoribo_table_core[i,]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid <- ath_idmapping$From[match(mitoribo_curr_row$Arabidopsis, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid)]
  
  mitoribo_curr_accessions <- c(mitoribo_curr_row$Mammals, mitoribo_curr_row$Yeast, mitoribo_curr_row$Trypanosomes, mitoribo_curr_row$Arabidopsis, mitoribo_curr_row$Spinacea_chloroplast_ribosome_6ERI, ath_uniprot_ids)
  
  # Only include homologs from structure
  OG_curr_row <- orthogroups_collapse_accessions_seprows_uniprot[which(orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession %in% mitoribo_curr_accessions),]
  OG_curr_row$mitoribo_id <- mitoribo_curr_row$NewName
  
  # Include homologs from other species
  OG_curr_all_members <- orthogroups_collapse_accessions_seprows_uniprot[orthogroups_collapse_accessions_seprows_uniprot$OG_id %in% OG_curr_row$OG_id,]
  giardia_index <- which(gsub("_.*", "", OG_curr_all_members$accessions) %in% c("5741"))
  if (length(giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[giardia_index,]
    OG_curr_all_members$mitoribo_id <- mitoribo_curr_row$NewName
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_mitoribo <- rbind(OG_mitoribo, OG_curr_row)
}

# Get taxids
OG_mitoribo$taxid <- gsub("_.*", "", OG_mitoribo$accessions)

phrogs_long_select_human <- phrogs_long %>% filter(protein_id %in% OG_mitoribo$accessions) %>% filter(taxid == "9606")
phrogs_long_select_mitoribo <- phrogs_long %>% filter(PROG_id %in% phrogs_long_select_human$PROG_id)

if (BOOL_USE_PHROG) {
  OG_mitoribo <- OG_mitoribo %>% filter(accessions %in% phrogs_long_select_mitoribo$protein_id)
} else {
  OG_mitoribo <- OG_mitoribo %>% group_by(mitoribo_id) %>% filter(OG_id %in% OG_id[taxid == "9606"]) %>% ungroup()
}

if (BOOL_USE_ABSENSE) {
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows, phrogs_long_select_human, by.x = c("human_accessions", "OG_id"), by.y = c("protein_id", "OG_id"))
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows_human, OG_mitoribo, by.x = c("human_accessions", "OG_id"), by.y = c("accessions", "OG_id"))
  homology_power_per_species_curr <- homology_power_per_species %>% filter(PhROG_id %in% homology_power_agg_seprows_human$PhROG_id) %>% filter(taxid %in% c("559292", "3702", "185431", "5741", "3562"))
  homology_power_per_species_curr$mitoribo_id <- homology_power_agg_seprows_human$mitoribo_id[match(homology_power_per_species_curr$PhROG_id, homology_power_agg_seprows_human$PhROG_id)]
  homology_power_per_species_human_underpowered <- homology_power_per_species_curr %>% filter(probability_of_detection < 0.95)
  OG_mitoribo_underpowered <- data.frame(OG_id = "Underpowered", accessions = NA, uniprot_accession = NA, mitoribo_id = homology_power_per_species_human_underpowered$mitoribo_id, taxid = homology_power_per_species_human_underpowered$taxid)
  OG_mitoribo <- rbind(OG_mitoribo, OG_mitoribo_underpowered)
  OG_mitoribo <- OG_mitoribo %>% group_by(mitoribo_id) %>% filter(!duplicated(taxid)) %>% ungroup()
}

OG_mitoribo_summary <- OG_mitoribo %>% group_by(mitoribo_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_mitoribo_species <- c()

for (mitoribo_id_curr in unique(OG_mitoribo_summary$mitoribo_id)) {
  OG_mitoribo_summary_curr <- OG_mitoribo_summary[OG_mitoribo_summary$mitoribo_id == mitoribo_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  giardia_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  spinach_chloroplast_ribosome_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_ribosome_OG_id) == 0) {
    spinach_chloroplast_ribosome_OG_id <- NA
  }
  
  OG_mitoribo_species_curr <- data.frame(mitoribo_id_curr, `9606` = human_OG_id, `559292` = yeast_OG_id, `185431` = trypanosoma_OG_id, `3702` = arabidopsis_OG_id, `5741` = giardia_OG_id, `3562` = spinach_chloroplast_ribosome_OG_id, check.names = FALSE)
  OG_mitoribo_species <- rbind(OG_mitoribo_species, OG_mitoribo_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Assign colors
if (BOOL_USE_ABSENSE) {
  OG_mitoribo_species_colors <- t(apply(OG_mitoribo_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", row == "Underpowered" ~ "gray", TRUE ~ "white")))
} else {
  OG_mitoribo_species_colors <- t(apply(OG_mitoribo_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", TRUE ~ "white")))
}
colnames(OG_mitoribo_species_colors) <- colnames(OG_mitoribo_species[,-1])
OG_mitoribo_species_colors <- as.data.frame(OG_mitoribo_species_colors)

# Reshape the color data into a dataframe for plotting
color_df <- OG_mitoribo_species_colors
color_df$row <- OG_mitoribo_species$mitoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(mitoribo_table_core$NewName))

if (BOOL_USE_PHROG) {
  dataset_description <- "PhROG"
} else {
  dataset_description <- "OG"
}

# Plotting the tile heatmap with ggplot2
pdf(paste0("mitoribo_", dataset_description, "_heatmap.pdf"), width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = dataset_description) +
  theme_minimal()
dev.off()




# Blast analysis
blast_out <- blast_out_BBH
title_string <- "BLAST BBH expect<1 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(blast_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
blast_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
blast_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(blast_out$sseqid, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

blast_out_mitoribo <- blast_out[which(blast_out$query_uniprot_accession %in% mitoribo_table_core$Mammals),]
blast_out_mitoribo$mitoribo_id <- mitoribo_table_core$NewName[match(blast_out_mitoribo$query_uniprot_accession, mitoribo_table_core$Mammals)]
blast_out_mitoribo <- blast_out_mitoribo %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% mitoribo_table_core[which(mitoribo_table_core$NewName == mitoribo_id), c("Mammals", "Yeast", "Trypanosomes", "Arabidopsis", "Spinacea_chloroplast_ribosome_6ERI")])
blast_out_mitoribo_summary <- blast_out_mitoribo %>% group_by(mitoribo_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Giardia hits as true, since don't have ground truth structure
blast_out_mitoribo_summary$any_true_top_hit[blast_out_mitoribo_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
blast_mitoribo_species <- c()
for (mitoribo_id_curr in unique(blast_out_mitoribo_summary$mitoribo_id)) {
  blast_out_mitoribo_curr <- blast_out_mitoribo_summary[blast_out_mitoribo_summary$mitoribo_id == mitoribo_id_curr,]
  
  # To use blast true hits
  blast_out_mitoribo_summary_curr_true_hit <- blast_out_mitoribo_curr[blast_out_mitoribo_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  blast_out_mitoribo_summary_curr <- data.frame(mitoribo_id_curr, human = 0, yeast = 0, trypanosoma = 0, arabidopsis = 0, giardia = 0, spinach_chloroplast = 0)
  
  # Binary
  blast_out_mitoribo_summary_curr$human <- as.numeric("9606" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid)
  blast_out_mitoribo_summary_curr$yeast <- as.numeric("559292" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid)
  blast_out_mitoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid)
  blast_out_mitoribo_summary_curr$arabidopsis <- as.numeric("3702" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid)
  blast_out_mitoribo_summary_curr$giardia <- as.numeric("5741" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid)
  blast_out_mitoribo_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid)
  
  # Tri color
  blast_out_mitoribo_summary_curr$human <- as.numeric("9606" %in% blast_out_mitoribo_curr$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$human <- blast_out_mitoribo_summary_curr$human + as.numeric("9606" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$yeast <- as.numeric("559292" %in% blast_out_mitoribo_curr$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$yeast <- blast_out_mitoribo_summary_curr$yeast + as.numeric("559292" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% blast_out_mitoribo_curr$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$trypanosoma <- blast_out_mitoribo_summary_curr$trypanosoma + as.numeric("185431" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$arabidopsis <- as.numeric("3702" %in% blast_out_mitoribo_curr$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$arabidopsis <- blast_out_mitoribo_summary_curr$arabidopsis + as.numeric("3702" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% blast_out_mitoribo_curr$target_taxid) * 0.5
  blast_out_mitoribo_summary_curr$spinach_chloroplast <- blast_out_mitoribo_summary_curr$spinach_chloroplast + as.numeric("3562" %in% blast_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  
  blast_mitoribo_species <- rbind(blast_mitoribo_species, blast_out_mitoribo_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(blast_mitoribo_species[,-1])

# Binary color
color_mat[color_mat == 0.5] <- 0
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)
color_df$row <- blast_mitoribo_species$mitoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(mitoribo_table_core$NewName))

# Plotting the tile heatmap with ggplot2
pdf("mitoribo_blast.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# plmsearch analysis
plmsearch_out <- plmsearch_out_BBH
title_string <- "plmsearch BBH similarity score > 0.3 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(plmsearch_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
plmsearch_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
plmsearch_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(plmsearch_out$target, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

plmsearch_out_mitoribo <- plmsearch_out[which(plmsearch_out$query_uniprot_accession %in% mitoribo_table_core$Mammals),]
plmsearch_out_mitoribo$mitoribo_id <- mitoribo_table_core$NewName[match(plmsearch_out_mitoribo$query_uniprot_accession, mitoribo_table_core$Mammals)]
plmsearch_out_mitoribo <- plmsearch_out_mitoribo %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% mitoribo_table_core[which(mitoribo_table_core$NewName == mitoribo_id), c("Mammals", "Yeast", "Trypanosomes", "Arabidopsis", "Spinacea_chloroplast_ribosome_6ERI")])
plmsearch_out_mitoribo_summary <- plmsearch_out_mitoribo %>% group_by(mitoribo_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Giardia hits as true, since don't have ground truth structure
plmsearch_out_mitoribo_summary$any_true_top_hit[plmsearch_out_mitoribo_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
plmsearch_mitoribo_species <- c()
for (mitoribo_id_curr in unique(plmsearch_out_mitoribo_summary$mitoribo_id)) {
  plmsearch_out_mitoribo_curr <- plmsearch_out_mitoribo_summary[plmsearch_out_mitoribo_summary$mitoribo_id == mitoribo_id_curr,]
  
  # To use plmsearch true hits
  plmsearch_out_mitoribo_summary_curr_true_hit <- plmsearch_out_mitoribo_curr[plmsearch_out_mitoribo_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  plmsearch_out_mitoribo_summary_curr <- data.frame(mitoribo_id_curr, human = 0, yeast = 0, trypanosoma = 0, arabidopsis = 0, giardia = 0, spinach_chloroplast = 0)
  
  # Binary
  plmsearch_out_mitoribo_summary_curr$human <- as.numeric("9606" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid)
  plmsearch_out_mitoribo_summary_curr$yeast <- as.numeric("559292" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid)
  plmsearch_out_mitoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid)
  plmsearch_out_mitoribo_summary_curr$arabidopsis <- as.numeric("3702" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid)
  plmsearch_out_mitoribo_summary_curr$giardia <- as.numeric("5741" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid)
  plmsearch_out_mitoribo_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid)
  
  # Tri color
  plmsearch_out_mitoribo_summary_curr$human <- as.numeric("9606" %in% plmsearch_out_mitoribo_curr$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$human <- plmsearch_out_mitoribo_summary_curr$human + as.numeric("9606" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$yeast <- as.numeric("559292" %in% plmsearch_out_mitoribo_curr$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$yeast <- plmsearch_out_mitoribo_summary_curr$yeast + as.numeric("559292" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% plmsearch_out_mitoribo_curr$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$trypanosoma <- plmsearch_out_mitoribo_summary_curr$trypanosoma + as.numeric("185431" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$arabidopsis <- as.numeric("3702" %in% plmsearch_out_mitoribo_curr$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$arabidopsis <- plmsearch_out_mitoribo_summary_curr$arabidopsis + as.numeric("3702" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% plmsearch_out_mitoribo_curr$target_taxid) * 0.5
  plmsearch_out_mitoribo_summary_curr$spinach_chloroplast <- plmsearch_out_mitoribo_summary_curr$spinach_chloroplast + as.numeric("3562" %in% plmsearch_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  
  plmsearch_mitoribo_species <- rbind(plmsearch_mitoribo_species, plmsearch_out_mitoribo_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(plmsearch_mitoribo_species[,-1])

# Binary color
color_mat[color_mat == 0.5] <- 0
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)
color_df$row <- plmsearch_mitoribo_species$mitoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(mitoribo_table_core$NewName))

# Plotting the tile heatmap with ggplot2
pdf("mitoribo_plmsearch.plmalign.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# Foldseek analysis
foldseek_out <- foldseek_out_BBH
title_string <- "foldseek BBH expect<1 with human homolog"

foldseek_out_mitoribo <- foldseek_out[which(foldseek_out$query %in% mitoribo_table_core$Mammals),]
foldseek_out_mitoribo$mitoribo_id <- mitoribo_table_core$NewName[match(foldseek_out_mitoribo$query, mitoribo_table_core$Mammals)]
# Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
foldseek_out_mitoribo <- foldseek_out_mitoribo %>% rowwise() %>% mutate(true_top_hit = target %in% c(mitoribo_table_core[which(mitoribo_table_core$NewName == mitoribo_id), c("Mammals", "Yeast", "Trypanosomes", "Arabidopsis", "Spinacea_chloroplast_ribosome_6ERI")], ath_idmapping$Entry[which(ath_idmapping$From == ath_idmapping$From[match(mitoribo_table_core$Arabidopsis[which(mitoribo_table_core$NewName == mitoribo_id)], ath_idmapping$Entry)])]))
foldseek_out_mitoribo_summary <- foldseek_out_mitoribo %>% group_by(mitoribo_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Giardia hits as true, since don't have ground truth structure
foldseek_out_mitoribo_summary$any_true_top_hit[foldseek_out_mitoribo_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
foldseek_mitoribo_species <- c()
for (mitoribo_id_curr in unique(foldseek_out_mitoribo_summary$mitoribo_id)) {
  foldseek_out_mitoribo_curr <- foldseek_out_mitoribo_summary[foldseek_out_mitoribo_summary$mitoribo_id == mitoribo_id_curr,]
  
  # To use foldseek true hits
  foldseek_out_mitoribo_summary_curr_true_hit <- foldseek_out_mitoribo_curr[foldseek_out_mitoribo_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  foldseek_out_mitoribo_summary_curr <- data.frame(mitoribo_id_curr, human = 0, yeast = 0, trypanosoma = 0, arabidopsis = 0, giardia = 0, spinach_chloroplast = 0)
  
  # Binary
  foldseek_out_mitoribo_summary_curr$human <- as.numeric("9606" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid)
  foldseek_out_mitoribo_summary_curr$yeast <- as.numeric("559292" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid)
  foldseek_out_mitoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid)
  foldseek_out_mitoribo_summary_curr$arabidopsis <- as.numeric("3702" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid)
  foldseek_out_mitoribo_summary_curr$giardia <- as.numeric("5741" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid)
  foldseek_out_mitoribo_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid)
  
  # Tri color
  foldseek_out_mitoribo_summary_curr$human <- as.numeric("9606" %in% foldseek_out_mitoribo_curr$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$human <- foldseek_out_mitoribo_summary_curr$human + as.numeric("9606" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$yeast <- as.numeric("559292" %in% foldseek_out_mitoribo_curr$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$yeast <- foldseek_out_mitoribo_summary_curr$yeast + as.numeric("559292" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% foldseek_out_mitoribo_curr$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$trypanosoma <- foldseek_out_mitoribo_summary_curr$trypanosoma + as.numeric("185431" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$arabidopsis <- as.numeric("3702" %in% foldseek_out_mitoribo_curr$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$arabidopsis <- foldseek_out_mitoribo_summary_curr$arabidopsis + as.numeric("3702" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% foldseek_out_mitoribo_curr$target_taxid) * 0.5
  foldseek_out_mitoribo_summary_curr$spinach_chloroplast <- foldseek_out_mitoribo_summary_curr$spinach_chloroplast + as.numeric("3562" %in% foldseek_out_mitoribo_summary_curr_true_hit$target_taxid) * 0.5
  
  foldseek_mitoribo_species <- rbind(foldseek_mitoribo_species, foldseek_out_mitoribo_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(foldseek_mitoribo_species[,-1])

# Binary color
color_mat[color_mat == 0.5] <- 0
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)
color_df$row <- foldseek_mitoribo_species$mitoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(mitoribo_table_core$NewName))

# Plotting the tile heatmap with ggplot2
pdf("mitoribo_foldseek.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




## EggNOG mitoribo
# Idmapping for human and yeast mitoribo
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "mitoribo_human_uniprot_to_ensemblprotein_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "mitoribo_yeast_uniprot_to_kegg_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

eggnog_species_ids <- c("9606", "4932", "3702", "185431")

OG_mitoribo <- c()
for (i in 1:nrow(mitoribo_table_core)) {
  mitoribo_curr_row <- mitoribo_table_core[i,]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid_curr <- ath_idmapping$From[match(mitoribo_curr_row$Arabidopsis, ath_idmapping$Entry)]
  ath_uniprot_ids <- paste0("3702.", ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid_curr)])
  
  # Map ids to eggnog
  ath_geneid <- paste0("3702.", mitoribo_curr_row$Arabidopsis)
  tbr_geneid <- paste0("185431.", mitoribo_curr_row$Trypanosomes)
  spinacea_geneid <- paste0("3562.", mitoribo_curr_row$Spinacea_chloroplast_ribosome_6ERI)
  sce_geneid <- paste0("4932.", idmapping_yeast$To[idmapping_yeast$From %in% mitoribo_curr_row$Yeast])
  hsa_geneid <- paste0("9606.", idmapping_human$To[idmapping_human$From %in% mitoribo_curr_row$Mammals])
  
  mitoribo_curr_accessions <- c(hsa_geneid, sce_geneid, ath_geneid, tbr_geneid, spinacea_geneid, ath_uniprot_ids)
  
  OG_curr_row <- root_OG_long %>% filter(protein_ids %in% mitoribo_curr_accessions)
  OG_curr_row$mitoribo_id <- mitoribo_curr_row$NewName
  
  # Include homologs from other species
  OG_curr_all_members <- root_OG_long[root_OG_long$OG_id %in% OG_curr_row$OG_id,]
  giardia_index <- which(gsub("\\..*", "", OG_curr_all_members$protein_ids) %in% c("184922"))
  if (length(giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[giardia_index,]
    OG_curr_all_members$mitoribo_id <- mitoribo_curr_row$NewName
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "protein_ids", "mitoribo_id")]
  
  # Find missing entries
  missing_species <- eggnog_species_ids[!eggnog_species_ids %in% gsub("\\..*", "", OG_curr_row$protein_ids)]
  if (length(missing_species) > 0) {
    print(mitoribo_curr_row$NewName)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", protein_ids = missing_species, mitoribo_id = mitoribo_curr_row$NewName)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_mitoribo <- rbind(OG_mitoribo, OG_curr_row)
}

# Get taxids
OG_mitoribo$taxid <- gsub("\\..*", "", OG_mitoribo$protein_ids)

OG_mitoribo_summary <- OG_mitoribo %>% group_by(mitoribo_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_mitoribo_species <- c()
for (mitoribo_id_curr in unique(OG_mitoribo_summary$mitoribo_id)) {
  OG_mitoribo_summary_curr <- OG_mitoribo_summary[OG_mitoribo_summary$mitoribo_id == mitoribo_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 4932]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  giardia_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 184922]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  spinach_chloroplast_ribosome_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_ribosome_OG_id) == 0) {
    spinach_chloroplast_ribosome_OG_id <- NA
  }
  
  OG_mitoribo_species_curr <- data.frame(mitoribo_id_curr, `9606` = human_OG_id, `559292` = yeast_OG_id, `185431` = trypanosoma_OG_id, `3702` = arabidopsis_OG_id, `5741` = giardia_OG_id, `3562` = spinach_chloroplast_ribosome_OG_id, check.names = FALSE)
  OG_mitoribo_species <- rbind(OG_mitoribo_species, OG_mitoribo_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_mitoribo_species_bool <- t(apply(OG_mitoribo_species[,-1], 1, function(row) row == row[1]))
OG_mitoribo_species_bool[is.na(OG_mitoribo_species_bool)] <- FALSE
OG_mitoribo_species_colors <- ifelse(OG_mitoribo_species_bool, "#0047AB", "white")
OG_mitoribo_species_colors <- as.data.frame(OG_mitoribo_species_colors)

# Mark proteins missing from dataset as gray
OG_mitoribo_species_colors[which(OG_mitoribo_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_mitoribo_species_colors

color_df$row <- OG_mitoribo_species$mitoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(mitoribo_table_core$NewName))


# Plotting the tile heatmap with ggplot2
pdf("mitoribo_EggNOG_Eukaryota_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "EggNOG") +
  theme_minimal()
dev.off()



## OrthoDB mitoribo
# Idmapping for human and yeast mitoribo
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "mitoribo_human_uniprot_to_ensembl_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "mitoribo_yeast_uniprot_to_kegg_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

orthodb_species_ids <- c("185431_0", "3702_0", "559292_0", "9606_0")

OG_mitoribo <- c()
for (i in 1:nrow(mitoribo_table_core)) {
  mitoribo_curr_row <- mitoribo_table_core[i,]
  
  sce_geneid <- idmapping_yeast$To[idmapping_yeast$From %in% mitoribo_curr_row$Yeast]
  hsa_ensg_id <- idmapping_human$To[idmapping_human$From %in% mitoribo_curr_row$Mammals]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid <- ath_idmapping$From[match(mitoribo_curr_row$Arabidopsis, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid)]
  idmapping_ath_curr <- idmapping_ath %>% filter(KEGG %in% ath_geneid | TAIR %in% ath_geneid | Araport %in% ath_geneid)
  ath_uniprot_ids <- c(ath_uniprot_ids, idmapping_ath_curr$From)
  
  tbr_gene_id <- idmapping_tbr$From[match(mitoribo_curr_row$Trypanosomes, idmapping_tbr$Entry)]
  
  mitoribo_curr_accessions <- c(mitoribo_curr_row$Trypanosomes, mitoribo_curr_row$Arabidopsis, mitoribo_curr_row$Spinacea_chloroplast_ribosome_6ERI, ath_uniprot_ids)
  mitoribo_curr_accessions <- mitoribo_curr_accessions[mitoribo_curr_accessions != ""]
  
  OG_curr_row <- odb_eukaryota_ogs_genes %>% filter(uniprot_id %in% mitoribo_curr_accessions | gene_id %in% hsa_ensg_id | synonym %in% sce_geneid | gene_id %in% ath_geneid | synonym %in% tbr_gene_id)
  OG_curr_row$mitoribo_id <- mitoribo_curr_row$NewName
  
  # Include homologs from other species
  OG_curr_all_members <- odb_eukaryota_ogs_genes[odb_eukaryota_ogs_genes$OG_id %in% OG_curr_row$OG_id,]
  giardia_index <- which(OG_curr_all_members$species_id == "5741_0")
  if (length(giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[giardia_index,]
    OG_curr_all_members$mitoribo_id <- mitoribo_curr_row$NewName
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "species_id", "mitoribo_id")]
  
  # Find missing entries
  missing_species <- orthodb_species_ids[!orthodb_species_ids %in% OG_curr_row$species_id]
  if (length(missing_species) > 0) {
    print(mitoribo_curr_row$NewName)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", species_id = missing_species, mitoribo_id = mitoribo_curr_row$NewName)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_mitoribo <- rbind(OG_mitoribo, OG_curr_row)
}

# Get taxids
OG_mitoribo$taxid <- gsub("_.*", "", OG_mitoribo$species_id)

OG_mitoribo_summary <- OG_mitoribo %>% group_by(mitoribo_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_mitoribo_species <- c()
for (mitoribo_id_curr in unique(OG_mitoribo_summary$mitoribo_id)) {
  OG_mitoribo_summary_curr <- OG_mitoribo_summary[OG_mitoribo_summary$mitoribo_id == mitoribo_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  giardia_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  spinach_chloroplast_ribosome_OG_id <- OG_mitoribo_summary_curr$OG_id[OG_mitoribo_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_ribosome_OG_id) == 0) {
    spinach_chloroplast_ribosome_OG_id <- NA
  }
  
  OG_mitoribo_species_curr <- data.frame(mitoribo_id_curr, `9606` = human_OG_id, `559292` = yeast_OG_id, `185431` = trypanosoma_OG_id, `3702` = arabidopsis_OG_id, `5741` = giardia_OG_id, `3562` = spinach_chloroplast_ribosome_OG_id, check.names = FALSE)
  OG_mitoribo_species <- rbind(OG_mitoribo_species, OG_mitoribo_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_mitoribo_species_bool <- t(apply(OG_mitoribo_species[,-1], 1, function(row) row == row[1]))
OG_mitoribo_species_bool[is.na(OG_mitoribo_species_bool)] <- FALSE
OG_mitoribo_species_colors <- ifelse(OG_mitoribo_species_bool, "#0047AB", "white")
OG_mitoribo_species_colors <- as.data.frame(OG_mitoribo_species_colors)

# Mark proteins missing from dataset as gray
OG_mitoribo_species_colors[which(OG_mitoribo_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_mitoribo_species_colors
color_df$row <- OG_mitoribo_species$mitoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(mitoribo_table_core$NewName))

# Plotting the tile heatmap with ggplot2
pdf("mitoribo_OrthoDB_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "OrthoDB") +
  theme_minimal()
dev.off()



## Complex V
complex5_table <- read.csv(here("data", "annotation", "pathways", "complexV_structure_annotation.csv"))

# Find conserved subunits. Exclude tetrahymena since it has some unmapped accessions (absent from uniprot reference proteome)
complex5_table <- complex5_table %>% rowwise() %>% mutate(n_spp = sum(c(!is.na(Homo_sapiens), !is.na(Saccharomyces_cerevisiae), !is.na(Trypanosoma_brucei), !is.na(Toxoplasma_gondii), !is.na(Arabidopsis_thaliana))))
complex5_table_core <- complex5_table[complex5_table$n_spp == 5,]

# OG data Complex V
complex_table_core <- complex5_table_core
complex_table_core$name <- complex_table_core$Subunit.Name..eukaryotic.nomenclature.
OG_complex <- c()
for (i in 1:nrow(complex_table_core)) {
  complex_curr_row <- complex_table_core[i,]
  
  complex_curr_accessions <- c(complex_curr_row$Homo_sapiens, complex_curr_row$Saccharomyces_cerevisiae, complex_curr_row$Trypanosoma_brucei, complex_curr_row$Tetrahymena_thermophila, complex_curr_row$Toxoplasma_gondii, complex_curr_row$Arabidopsis_thaliana, complex_curr_row$Spinacia_olacrea_chloroplast, complex_curr_row$Escherichia_coli)
  
  # Only include homologs from structure
  OG_curr_row <- orthogroups_collapse_accessions_seprows_uniprot[which(orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession %in% complex_curr_accessions),]
  
  # Include homologs from Giardia
  OG_curr_all_members <- orthogroups_collapse_accessions_seprows_uniprot[orthogroups_collapse_accessions_seprows_uniprot$OG_id %in% OG_curr_row$OG_id,]
  giardia_index <- which(gsub("_.*", "", OG_curr_all_members$accessions) %in% c("5741"))
  if (length(giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[giardia_index,]
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row$complex_id <- complex_curr_row$name
  
  OG_complex <- rbind(OG_complex, OG_curr_row)
}

# Get taxids
OG_complex$taxid <- gsub("_.*", "", OG_complex$accessions)

phrogs_long_select_human <- phrogs_long %>% filter(protein_id %in% OG_complex$accessions) %>% filter(taxid == "9606")
phrogs_long_select_complex <- phrogs_long %>% filter(PROG_id %in% phrogs_long_select_human$PROG_id)

if (BOOL_USE_PHROG) {
  OG_complex <- OG_complex %>% filter(accessions %in% phrogs_long_select_complex$protein_id)
} else {
  OG_complex <- OG_complex %>% group_by(complex_id) %>% filter(OG_id %in% OG_id[taxid == "9606"]) %>% ungroup()
}

if (BOOL_USE_ABSENSE) {
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows, phrogs_long_select_human, by.x = c("human_accessions", "OG_id"), by.y = c("protein_id", "OG_id"))
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows_human, OG_complex, by.x = c("human_accessions", "OG_id"), by.y = c("accessions", "OG_id"))
  homology_power_per_species_curr <- homology_power_per_species %>% filter(PhROG_id %in% homology_power_agg_seprows_human$PhROG_id) %>% filter(taxid %in% c("559292", "3702", "185431", "508771", "5741", "3562"))
  homology_power_per_species_curr$complex_id <- homology_power_agg_seprows_human$complex_id[match(homology_power_per_species_curr$PhROG_id, homology_power_agg_seprows_human$PhROG_id)]
  homology_power_per_species_human_underpowered <- homology_power_per_species_curr %>% filter(probability_of_detection < 0.95)
  OG_complex_underpowered <- data.frame(OG_id = "Underpowered", accessions = NA, uniprot_accession = NA, complex_id = homology_power_per_species_human_underpowered$complex_id, taxid = homology_power_per_species_human_underpowered$taxid)
  OG_complex <- rbind(OG_complex, OG_complex_underpowered)
  OG_complex <- OG_complex %>% group_by(complex_id) %>% filter(!duplicated(taxid)) %>% ungroup()
}

OG_complex_summary <- OG_complex %>% group_by(complex_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_complex_species <- c()
for (complex_id_curr in unique(OG_complex_summary$complex_id)) {
  OG_complex_summary_curr <- OG_complex_summary[OG_complex_summary$complex_id == complex_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  toxoplasma_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 508771]
  if (length(toxoplasma_OG_id) == 0) {
    toxoplasma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  spinach_chloroplast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_OG_id) == 0) {
    spinach_chloroplast_OG_id <- NA
  }
  giardia_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  
  OG_complex_species_curr <- data.frame(complex_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Toxoplasma = toxoplasma_OG_id, Arabidopsis = arabidopsis_OG_id, Giardia = giardia_OG_id, Spinach_chloroplast = spinach_chloroplast_OG_id)
  OG_complex_species <- rbind(OG_complex_species, OG_complex_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen", "brown", "black") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Assign colors
if (BOOL_USE_ABSENSE) {
  OG_complex_species_colors <- t(apply(OG_complex_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", row == "Underpowered" ~ "gray", TRUE ~ "white")))
} else {
  OG_complex_species_colors <- t(apply(OG_complex_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", TRUE ~ "white")))
}
colnames(OG_complex_species_colors) <- colnames(OG_complex_species[,-1])
OG_complex_species_colors <- as.data.frame(OG_complex_species_colors)

# Reshape the color data into a dataframe for plotting
color_df <- OG_complex_species_colors

color_df$row <- OG_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$name))

if (BOOL_USE_PHROG) {
  dataset_description <- "PhROG"
} else {
  dataset_description <- "OG"
}

# Plotting the tile heatmap with ggplot2
pdf(paste0("complexV_", dataset_description, "_heatmap.pdf"), width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = dataset_description) +
  theme_minimal()
dev.off()




# Blast analysis
blast_out <- blast_out_BBH
title_string <- "BLAST BBH expect<1 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(blast_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
blast_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
blast_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(blast_out$sseqid, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

blast_out_complex5 <- blast_out[which(blast_out$query_uniprot_accession %in% complex5_table_core$Homo_sapiens),]
blast_out_complex5$complex5_id <- complex5_table_core$Subunit.Name..eukaryotic.nomenclature.[match(blast_out_complex5$query_uniprot_accession, complex5_table_core$Homo_sapiens)]
blast_out_complex5 <- blast_out_complex5 %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% complex5_table_core[which(complex5_table_core$Subunit.Name..eukaryotic.nomenclature. == complex5_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Toxoplasma_gondii", "Arabidopsis_thaliana", "Spinacia_olacrea_chloroplast")])
blast_out_complex5_summary <- blast_out_complex5 %>% group_by(complex5_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Giardia hits as true, since don't have ground truth structure
blast_out_complex5_summary$any_true_top_hit[blast_out_complex5_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
blast_complex5_species <- c()
for (complex5_id_curr in unique(blast_out_complex5_summary$complex5_id)) {
  blast_out_complex5_curr <- blast_out_complex5_summary[blast_out_complex5_summary$complex5_id == complex5_id_curr,]
  
  # To use blast true hits
  blast_out_complex5_summary_curr_true_hit <- blast_out_complex5_curr[blast_out_complex5_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  blast_out_complex5_summary_curr <- data.frame(complex5_id_curr, Human = 0, Yeast = 0, Trypanosoma = 0, Toxoplasma = 0, Arabidopsis = 0, Giardia = 0, Spinach_chloroplast = 0)
  
  # Binary
  blast_out_complex5_summary_curr$Human <- as.numeric("9606" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  blast_out_complex5_summary_curr$Yeast <- as.numeric("559292" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  blast_out_complex5_summary_curr$Trypanosoma <- as.numeric("185431" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  blast_out_complex5_summary_curr$Toxoplasma <- as.numeric("508771" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  blast_out_complex5_summary_curr$Arabidopsis <- as.numeric("3702" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  blast_out_complex5_summary_curr$Giardia <- as.numeric("5741" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  blast_out_complex5_summary_curr$Spinach_chloroplast <- as.numeric("3562" %in% blast_out_complex5_summary_curr_true_hit$target_taxid)
  
  blast_complex5_species <- rbind(blast_complex5_species, blast_out_complex5_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(blast_complex5_species[,-1])

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- blast_complex5_species$complex5_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex5_table_core$Subunit.Name..eukaryotic.nomenclature.))

# Plotting the tile heatmap with ggplot2
pdf("complexV_blast.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# plmsearch analysis
plmsearch_out <- plmsearch_out_BBH
title_string <- "plmsearch BBH similarity score > 0.3 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(plmsearch_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
plmsearch_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
plmsearch_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(plmsearch_out$target, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

plmsearch_out_complex5 <- plmsearch_out[which(plmsearch_out$query_uniprot_accession %in% complex5_table_core$Homo_sapiens),]
plmsearch_out_complex5$complex5_id <- complex5_table_core$Subunit.Name..eukaryotic.nomenclature.[match(plmsearch_out_complex5$query_uniprot_accession, complex5_table_core$Homo_sapiens)]
plmsearch_out_complex5 <- plmsearch_out_complex5 %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% complex5_table_core[which(complex5_table_core$Subunit.Name..eukaryotic.nomenclature. == complex5_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Toxoplasma_gondii", "Arabidopsis_thaliana", "Spinacia_olacrea_chloroplast")])
plmsearch_out_complex5_summary <- plmsearch_out_complex5 %>% group_by(complex5_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Giardia hits as true, since don't have ground truth structure
plmsearch_out_complex5_summary$any_true_top_hit[plmsearch_out_complex5_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
plmsearch_complex5_species <- c()
for (complex5_id_curr in unique(plmsearch_out_complex5_summary$complex5_id)) {
  plmsearch_out_complex5_curr <- plmsearch_out_complex5_summary[plmsearch_out_complex5_summary$complex5_id == complex5_id_curr,]
  
  # To use plmsearch true hits
  plmsearch_out_complex5_summary_curr_true_hit <- plmsearch_out_complex5_curr[plmsearch_out_complex5_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  plmsearch_out_complex5_summary_curr <- data.frame(complex5_id_curr, Human = 0, Yeast = 0, Trypanosoma = 0, Toxoplasma = 0, Arabidopsis = 0, Giardia = 0, Spinach_chloroplast = 0)
  
  # Binary
  plmsearch_out_complex5_summary_curr$Human <- as.numeric("9606" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex5_summary_curr$Yeast <- as.numeric("559292" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex5_summary_curr$Trypanosoma <- as.numeric("185431" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex5_summary_curr$Toxoplasma <- as.numeric("508771" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex5_summary_curr$Arabidopsis <- as.numeric("3702" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex5_summary_curr$Giardia <- as.numeric("5741" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex5_summary_curr$Spinach_chloroplast <- as.numeric("3562" %in% plmsearch_out_complex5_summary_curr_true_hit$target_taxid)
  
  plmsearch_complex5_species <- rbind(plmsearch_complex5_species, plmsearch_out_complex5_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(plmsearch_complex5_species[,-1])

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")
color_df <- as.data.frame(color_mat)

color_df$row <- plmsearch_complex5_species$complex5_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex5_table_core$Subunit.Name..eukaryotic.nomenclature.))

# Plotting the tile heatmap with ggplot2
pdf("complexV_plmsearch.plmalign.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# foldseek analysis
foldseek_out <- foldseek_out_BBH
title_string <- "foldseek BBH expect<1 with human homolog"

foldseek_out_complex5 <- foldseek_out[which(foldseek_out$query %in% complex5_table_core$Homo_sapiens),]
foldseek_out_complex5$complex5_id <- complex5_table_core$Subunit.Name..eukaryotic.nomenclature.[match(foldseek_out_complex5$query, complex5_table_core$Homo_sapiens)]
# Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
foldseek_out_complex5 <- foldseek_out_complex5 %>% rowwise() %>% mutate(true_top_hit = target %in% c(complex5_table_core[which(complex5_table_core$Subunit.Name..eukaryotic.nomenclature. == complex5_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Toxoplasma_gondii", "Arabidopsis_thaliana", "Spinacia_olacrea_chloroplast")], ath_idmapping$Entry[which(ath_idmapping$From == ath_idmapping$From[match(complex5_table_core$Arabidopsis_thaliana[which(complex5_table_core$Subunit.Name..eukaryotic.nomenclature. == complex5_id)], ath_idmapping$Entry)])]))
foldseek_out_complex5_summary <- foldseek_out_complex5 %>% group_by(complex5_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Giardia hits as true, since don't have ground truth structure
foldseek_out_complex5_summary$any_true_top_hit[foldseek_out_complex5_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
foldseek_complex5_species <- c()
for (complex5_id_curr in unique(foldseek_out_complex5_summary$complex5_id)) {
  foldseek_out_complex5_curr <- foldseek_out_complex5_summary[foldseek_out_complex5_summary$complex5_id == complex5_id_curr,]
  
  # To use foldseek true hits
  foldseek_out_complex5_summary_curr_true_hit <- foldseek_out_complex5_curr[foldseek_out_complex5_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  foldseek_out_complex5_summary_curr <- data.frame(complex5_id_curr, Human = 0, Yeast = 0, Trypanosoma = 0, Toxoplasma = 0, Arabidopsis = 0, Giardia = 0, Spinach_chloroplast = 0)
  
  # Binary
  foldseek_out_complex5_summary_curr$Human <- as.numeric("9606" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  foldseek_out_complex5_summary_curr$Yeast <- as.numeric("559292" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  foldseek_out_complex5_summary_curr$Trypanosoma <- as.numeric("185431" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  foldseek_out_complex5_summary_curr$Toxoplasma <- as.numeric("508771" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  foldseek_out_complex5_summary_curr$Arabidopsis <- as.numeric("3702" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  foldseek_out_complex5_summary_curr$Giardia <- as.numeric("5741" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  foldseek_out_complex5_summary_curr$Spinach_chloroplast <- as.numeric("3562" %in% foldseek_out_complex5_summary_curr_true_hit$target_taxid)
  
  foldseek_complex5_species <- rbind(foldseek_complex5_species, foldseek_out_complex5_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(foldseek_complex5_species[,-1])

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- foldseek_complex5_species$complex5_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex5_table_core$Subunit.Name..eukaryotic.nomenclature.))

# Plotting the tile heatmap with ggplot2
pdf("complexV_foldseek.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




## EggNOG complex V
# Idmapping for human and yeast complex V
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "complexV_human_uniprot_to_ensemblprotein_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "complexV_yeast_uniprot_to_kegg_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

eggnog_species_ids <- c("9606", "4932", "3702", "185431", "508771", "3562")

OG_complex <- c()
for (i in 1:nrow(complex_table_core)) {
  complex_curr_row <- complex_table_core[i,]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid_curr <- ath_idmapping$From[match(complex_curr_row$Arabidopsis_thaliana, ath_idmapping$Entry)]
  ath_uniprot_ids <- paste0("3702.", ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid_curr)])
  
  # Map ids to eggnog
  ath_geneid <- paste0("3702.", complex_curr_row$Arabidopsis_thaliana)
  tbr_geneid <- paste0("185431.", complex_curr_row$Trypanosoma_brucei)
  toxo_geneid <- paste0("508771.", complex_curr_row$Toxoplasma_gondii)
  sce_geneid <- paste0("4932.", idmapping_yeast$To[idmapping_yeast$From %in% complex_curr_row$Saccharomyces_cerevisiae])
  hsa_geneid <- paste0("9606.", idmapping_human$To[idmapping_human$From %in% complex_curr_row$Homo_sapiens])
  spinach_chloroplast_geneid <- paste0("3562.", complex_curr_row$Spinacia_olacrea_chloroplast)
  
  complex_curr_accessions <- c(hsa_geneid, sce_geneid, ath_geneid, tbr_geneid, toxo_geneid, spinach_chloroplast_geneid, ath_uniprot_ids)
  
  OG_curr_row <- root_OG_long %>% filter(protein_ids %in% complex_curr_accessions)
  OG_curr_row$complex_id <- complex_curr_row$name
  
  # Include homologs from other species
  OG_curr_all_members <- root_OG_long[root_OG_long$OG_id %in% OG_curr_row$OG_id,]
  giardia_index <- which(gsub("\\..*", "", OG_curr_all_members$protein_ids) %in% c("184922"))
  if (length(giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[giardia_index,]
    OG_curr_all_members$complex_id <- complex_curr_row$name
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "protein_ids", "complex_id")]
  
  # Find missing entries
  missing_species <- eggnog_species_ids[!eggnog_species_ids %in% gsub("\\..*", "", OG_curr_row$protein_ids)]
  if (length(missing_species) > 0) {
    print(complex_curr_row$name)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", protein_ids = missing_species, complex_id = complex_curr_row$name)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_complex <- rbind(OG_complex, OG_curr_row)
}

# Get taxids
OG_complex$taxid <- gsub("\\..*", "", OG_complex$protein_ids)

OG_complex_summary <- OG_complex %>% group_by(complex_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_complex_species <- c()
for (complex_id_curr in unique(OG_complex_summary$complex_id)) {
  OG_complex_summary_curr <- OG_complex_summary[OG_complex_summary$complex_id == complex_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 4932]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  toxoplasma_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 508771]
  if (length(toxoplasma_OG_id) == 0) {
    toxoplasma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  spinach_chloroplast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_OG_id) == 0) {
    spinach_chloroplast_OG_id <- NA
  }
  giardia_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 184922]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }

  OG_complex_species_curr <- data.frame(complex_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Toxoplasma = toxoplasma_OG_id, Arabidopsis = arabidopsis_OG_id, Giardia = giardia_OG_id, Spinach_chloroplast = spinach_chloroplast_OG_id)
  OG_complex_species <- rbind(OG_complex_species, OG_complex_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_complex_species_bool <- t(apply(OG_complex_species[,-1], 1, function(row) row == row[1]))
OG_complex_species_bool[is.na(OG_complex_species_bool)] <- FALSE
OG_complex_species_colors <- ifelse(OG_complex_species_bool, "#0047AB", "white")
OG_complex_species_colors <- as.data.frame(OG_complex_species_colors)

# Mark proteins missing from dataset as gray
OG_complex_species_colors[which(OG_complex_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_complex_species_colors

color_df$row <- OG_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$name))


# Plotting the tile heatmap with ggplot2
pdf("complexV_EggNOG_Eukaryota_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "EggNOG") +
  theme_minimal()
dev.off()



## OrthoDB Complex V

# Idmapping for human and yeast Complex V
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "complexV_human_uniprot_to_ensembl_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "complexV_yeast_uniprot_to_kegg_idmapping_2025_07_02.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

orthodb_species_ids <- c("185431_0", "3702_0", "559292_0", "9606_0", "508771_0", "3562_0")

complex_table_core <- complex5_table_core
complex_table_core$name <- complex_table_core$Subunit.Name..eukaryotic.nomenclature.
OG_complex <- c()
for (i in 1:nrow(complex_table_core)) {
  complex_curr_row <- complex_table_core[i,]
  
  sce_geneid <- idmapping_yeast$To[idmapping_yeast$From %in% complex_curr_row$Saccharomyces_cerevisiae]
  hsa_ensg_id <- idmapping_human$To[idmapping_human$From %in% complex_curr_row$Homo_sapiens]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid <- ath_idmapping$From[match(complex_curr_row$Arabidopsis_thaliana, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid)]
  idmapping_ath_curr <- idmapping_ath %>% filter(KEGG %in% ath_geneid | TAIR %in% ath_geneid | Araport %in% ath_geneid)
  ath_uniprot_ids <- c(ath_uniprot_ids, idmapping_ath_curr$From)
  
  tbr_gene_id <- idmapping_tbr$From[match(complex_curr_row$Trypanosoma_brucei, idmapping_tbr$Entry)]
  
  complex_curr_accessions <- c(complex_curr_row$Trypanosoma_brucei, complex_curr_row$Arabidopsis_thaliana, complex_curr_row$Toxoplasma_gondii, complex_curr_row$Spinacia_olacrea_chloroplast, ath_uniprot_ids)
  complex_curr_accessions <- complex_curr_accessions[complex_curr_accessions != ""]
  
  OG_curr_row <- odb_eukaryota_ogs_genes %>% filter(uniprot_id %in% complex_curr_accessions | gene_id %in% hsa_ensg_id | synonym %in% sce_geneid | gene_id %in% ath_geneid | synonym %in% tbr_gene_id)
  OG_curr_row$complex_id <- complex_curr_row$name
  
  # Include homologs from other species
  OG_curr_all_members <- odb_eukaryota_ogs_genes[odb_eukaryota_ogs_genes$OG_id %in% OG_curr_row$OG_id,]
  giardia_index <- which(OG_curr_all_members$species_id == "5741_0")
  if (length(giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[giardia_index,]
    OG_curr_all_members$complex_id <- complex_curr_row$name
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "species_id", "complex_id")]
  
  # Find missing entries
  missing_species <- orthodb_species_ids[!orthodb_species_ids %in% OG_curr_row$species_id]
  if (length(missing_species) > 0) {
    print(complex_curr_row$name)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", species_id = missing_species, complex_id = complex_curr_row$name)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_complex <- rbind(OG_complex, OG_curr_row)
}

# Get taxids
OG_complex$taxid <- gsub("_.*", "", OG_complex$species_id)

OG_complex_summary <- OG_complex %>% group_by(complex_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_complex_species <- c()
for (complex_id_curr in unique(OG_complex_summary$complex_id)) {
  OG_complex_summary_curr <- OG_complex_summary[OG_complex_summary$complex_id == complex_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  toxoplasma_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 508771]
  if (length(toxoplasma_OG_id) == 0) {
    toxoplasma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  spinach_chloroplast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_OG_id) == 0) {
    spinach_chloroplast_OG_id <- NA
  }
  giardia_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  
  OG_complex_species_curr <- data.frame(complex_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Toxoplasma = toxoplasma_OG_id, Arabidopsis = arabidopsis_OG_id, Giardia = giardia_OG_id, Spinach_chloroplast = spinach_chloroplast_OG_id)
  OG_complex_species <- rbind(OG_complex_species, OG_complex_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_complex_species_bool <- t(apply(OG_complex_species[,-1], 1, function(row) row == row[1]))
OG_complex_species_bool[is.na(OG_complex_species_bool)] <- FALSE
OG_complex_species_colors <- ifelse(OG_complex_species_bool, "#0047AB", "white")
OG_complex_species_colors <- as.data.frame(OG_complex_species_colors)

# Mark proteins missing from dataset as gray
OG_complex_species_colors[which(OG_complex_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_complex_species_colors

color_df$row <- OG_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$name))


# Plotting the tile heatmap with ggplot2
pdf("complexV_OrthoDB_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "OrthoDB") +
  theme_minimal()
dev.off()


## Cytoribo
# Read in cytoribo gene sets
cytoribo_table <- read.csv(here("data", "annotation", "pathways", "cytoribo_structure_annotation.csv"), header=TRUE)

cytoribo_table_core <- cytoribo_table %>% rowwise() %>% filter(!is.na(Homo_sapiens) & !is.na(Saccharomyces_cerevisiae) & !is.na(Trypanosoma_brucei) & !is.na(Giardia_intestinalis))


# OG data cytoribo
OG_cytoribo <- c()
for (i in 1:nrow(cytoribo_table_core)) {
  cytoribo_curr_row <- cytoribo_table_core[i,]
  
  cytoribo_curr_accessions <- c(cytoribo_curr_row$Homo_sapiens, cytoribo_curr_row$Saccharomyces_cerevisiae, cytoribo_curr_row$Trypanosoma_brucei, cytoribo_curr_row$Giardia_intestinalis)
  
  # Only include homologs from structure
  OG_curr_row <- orthogroups_collapse_accessions_seprows_uniprot[which(orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession %in% cytoribo_curr_accessions),]
  OG_curr_row$cytoribo_id <- cytoribo_curr_row$Subunit.Name
  
  OG_cytoribo <- rbind(OG_cytoribo, OG_curr_row)
}

# Get taxids
OG_cytoribo$taxid <- gsub("_.*", "", OG_cytoribo$accessions)

phrogs_long_select_human <- phrogs_long %>% filter(protein_id %in% OG_cytoribo$accessions) %>% filter(taxid == "9606")
phrogs_long_select_cytoribo <- phrogs_long %>% filter(PROG_id %in% phrogs_long_select_human$PROG_id)

if (BOOL_USE_PHROG) {
  OG_cytoribo <- OG_cytoribo %>% filter(accessions %in% phrogs_long_select_cytoribo$protein_id)
} else {
  OG_cytoribo <- OG_cytoribo %>% group_by(cytoribo_id) %>% filter(OG_id %in% OG_id[taxid == "9606"]) %>% ungroup()
}

if (BOOL_USE_ABSENSE) {
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows, phrogs_long_select_human, by.x = c("human_accessions", "OG_id"), by.y = c("protein_id", "OG_id"))
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows_human, OG_cytoribo, by.x = c("human_accessions", "OG_id"), by.y = c("accessions", "OG_id"))
  homology_power_per_species_curr <- homology_power_per_species %>% filter(PhROG_id %in% homology_power_agg_seprows_human$PhROG_id) %>% filter(taxid %in% c("559292", "185431", "5741"))
  homology_power_per_species_curr$cytoribo_id <- homology_power_agg_seprows_human$cytoribo_id[match(homology_power_per_species_curr$PhROG_id, homology_power_agg_seprows_human$PhROG_id)]
  homology_power_per_species_human_underpowered <- homology_power_per_species_curr %>% filter(probability_of_detection < 0.95)
  if (nrow(homology_power_per_species_human_underpowered) > 0) {
    OG_cytoribo_underpowered <- data.frame(OG_id = "Underpowered", accessions = NA, uniprot_accession = NA, cytoribo_id = homology_power_per_species_human_underpowered$cytoribo_id, taxid = homology_power_per_species_human_underpowered$taxid)
    OG_cytoribo <- rbind(OG_cytoribo, OG_cytoribo_underpowered)
    OG_cytoribo <- OG_cytoribo %>% group_by(cytoribo_id) %>% filter(!duplicated(taxid)) %>% ungroup()
  }
}

OG_cytoribo_summary <- OG_cytoribo %>% group_by(cytoribo_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_cytoribo_species <- c()
for (cytoribo_id_curr in unique(OG_cytoribo_summary$cytoribo_id)) {
  OG_cytoribo_summary_curr <- OG_cytoribo_summary[OG_cytoribo_summary$cytoribo_id == cytoribo_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  giardia_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  
  OG_cytoribo_species_curr <- data.frame(cytoribo_id_curr, human = human_OG_id, yeast = yeast_OG_id, trypanosoma = trypanosoma_OG_id, giardia = giardia_OG_id)
  OG_cytoribo_species <- rbind(OG_cytoribo_species, OG_cytoribo_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Assign colors
if (BOOL_USE_ABSENSE) {
  OG_cytoribo_species_colors <- t(apply(OG_cytoribo_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", row == "Underpowered" ~ "gray", TRUE ~ "white")))
} else {
  OG_cytoribo_species_colors <- t(apply(OG_cytoribo_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", TRUE ~ "white")))
}
colnames(OG_cytoribo_species_colors) <- colnames(OG_cytoribo_species[,-1])
OG_cytoribo_species_colors <- as.data.frame(OG_cytoribo_species_colors)

# Reshape the color data into a dataframe for plotting
color_df <- OG_cytoribo_species_colors
color_df$row <- OG_cytoribo_species$cytoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(cytoribo_table_core$Subunit.Name))

if (BOOL_USE_PHROG) {
  dataset_description <- "PhROG"
} else {
  dataset_description <- "OG"
}

# Plotting the tile heatmap with ggplot2
pdf(paste0("cytoribo_", dataset_description, "_heatmap.pdf"), width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = dataset_description) +
  theme_minimal()
dev.off()




# Blast analysis
# Read in blast BBH
blast_out <- blast_out_BBH
title_string <- "BLAST BBH expect<1 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(blast_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
blast_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
blast_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(blast_out$sseqid, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

blast_out_cytoribo <- blast_out[which(blast_out$query_uniprot_accession %in% cytoribo_table_core$Homo_sapiens),]
blast_out_cytoribo$cytoribo_id <- cytoribo_table_core$Subunit.Name[match(blast_out_cytoribo$query_uniprot_accession, cytoribo_table_core$Homo_sapiens)]
blast_out_cytoribo <- blast_out_cytoribo %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% cytoribo_table_core[which(cytoribo_table_core$Subunit.Name == cytoribo_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Giardia_intestinalis")])
blast_out_cytoribo_summary <- blast_out_cytoribo %>% group_by(cytoribo_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Split into 1 column per species
blast_cytoribo_species <- c()
for (cytoribo_id_curr in unique(blast_out_cytoribo_summary$cytoribo_id)) {
  blast_out_cytoribo_curr <- blast_out_cytoribo_summary[blast_out_cytoribo_summary$cytoribo_id == cytoribo_id_curr,]
  
  # To use blast true hits
  blast_out_cytoribo_summary_curr_true_hit <- blast_out_cytoribo_curr[blast_out_cytoribo_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  blast_out_cytoribo_summary_curr <- data.frame(cytoribo_id_curr, human = 0, yeast = 0, trypanosoma = 0, giardia = 0)
  blast_out_cytoribo_summary_curr$human <- as.numeric("9606" %in% blast_out_cytoribo_curr$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$human <- blast_out_cytoribo_summary_curr$human + as.numeric("9606" %in% blast_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$yeast <- as.numeric("559292" %in% blast_out_cytoribo_curr$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$yeast <- blast_out_cytoribo_summary_curr$yeast + as.numeric("559292" %in% blast_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% blast_out_cytoribo_curr$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$trypanosoma <- blast_out_cytoribo_summary_curr$trypanosoma + as.numeric("185431" %in% blast_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$giardia <- as.numeric("5741" %in% blast_out_cytoribo_curr$target_taxid) * 0.5
  blast_out_cytoribo_summary_curr$giardia <- blast_out_cytoribo_summary_curr$giardia + as.numeric("5741" %in% blast_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  
  blast_cytoribo_species <- rbind(blast_cytoribo_species, blast_out_cytoribo_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(blast_cytoribo_species[,-1])
color_mat[color_mat == 0.5] <- 0

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- blast_cytoribo_species$cytoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(cytoribo_table_core$Subunit.Name))

# Plotting the tile heatmap with ggplot2
pdf("cytoribo_blast.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# plmsearch analysis
# Read in plmsearch BBH
plmsearch_out <- plmsearch_out_BBH
title_string <- "plmsearch BBH similarity score > 0.3 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(plmsearch_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
plmsearch_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
plmsearch_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(plmsearch_out$target, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

plmsearch_out_cytoribo <- plmsearch_out[which(plmsearch_out$query_uniprot_accession %in% cytoribo_table_core$Homo_sapiens),]
plmsearch_out_cytoribo$cytoribo_id <- cytoribo_table_core$Subunit.Name[match(plmsearch_out_cytoribo$query_uniprot_accession, cytoribo_table_core$Homo_sapiens)]
plmsearch_out_cytoribo <- plmsearch_out_cytoribo %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% cytoribo_table_core[which(cytoribo_table_core$Subunit.Name == cytoribo_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Giardia_intestinalis")])
plmsearch_out_cytoribo_summary <- plmsearch_out_cytoribo %>% group_by(cytoribo_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Split into 1 column per species
plmsearch_cytoribo_species <- c()
for (cytoribo_id_curr in unique(plmsearch_out_cytoribo_summary$cytoribo_id)) {
  plmsearch_out_cytoribo_curr <- plmsearch_out_cytoribo_summary[plmsearch_out_cytoribo_summary$cytoribo_id == cytoribo_id_curr,]
  
  # To use plmsearch true hits
  plmsearch_out_cytoribo_summary_curr_true_hit <- plmsearch_out_cytoribo_curr[plmsearch_out_cytoribo_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  plmsearch_out_cytoribo_summary_curr <- data.frame(cytoribo_id_curr, human = 0, yeast = 0, trypanosoma = 0, giardia = 0)
  plmsearch_out_cytoribo_summary_curr$human <- as.numeric("9606" %in% plmsearch_out_cytoribo_curr$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$human <- plmsearch_out_cytoribo_summary_curr$human + as.numeric("9606" %in% plmsearch_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$yeast <- as.numeric("559292" %in% plmsearch_out_cytoribo_curr$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$yeast <- plmsearch_out_cytoribo_summary_curr$yeast + as.numeric("559292" %in% plmsearch_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% plmsearch_out_cytoribo_curr$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$trypanosoma <- plmsearch_out_cytoribo_summary_curr$trypanosoma + as.numeric("185431" %in% plmsearch_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$giardia <- as.numeric("5741" %in% plmsearch_out_cytoribo_curr$target_taxid) * 0.5
  plmsearch_out_cytoribo_summary_curr$giardia <- plmsearch_out_cytoribo_summary_curr$giardia + as.numeric("5741" %in% plmsearch_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  
  plmsearch_cytoribo_species <- rbind(plmsearch_cytoribo_species, plmsearch_out_cytoribo_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(plmsearch_cytoribo_species[,-1])
color_mat[color_mat == 0.5] <- 0

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- plmsearch_cytoribo_species$cytoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(cytoribo_table_core$Subunit.Name))

# Plotting the tile heatmap with ggplot2
pdf("cytoribo_plmsearch.plmalign.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# foldseek analysis
foldseek_out <- foldseek_out_BBH
title_string <- "foldseek BBH expect<1 with human homolog"

foldseek_out_cytoribo <- foldseek_out[which(foldseek_out$query %in% cytoribo_table_core$Homo_sapiens),]
foldseek_out_cytoribo$cytoribo_id <- cytoribo_table_core$Subunit.Name[match(foldseek_out_cytoribo$query, cytoribo_table_core$Homo_sapiens)]
foldseek_out_cytoribo <- foldseek_out_cytoribo %>% rowwise() %>% mutate(true_top_hit = target %in% cytoribo_table_core[which(cytoribo_table_core$Subunit.Name == cytoribo_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Giardia_intestinalis")])
foldseek_out_cytoribo_summary <- foldseek_out_cytoribo %>% group_by(cytoribo_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Split into 1 column per species
foldseek_cytoribo_species <- c()
for (cytoribo_id_curr in unique(foldseek_out_cytoribo_summary$cytoribo_id)) {
  foldseek_out_cytoribo_curr <- foldseek_out_cytoribo_summary[foldseek_out_cytoribo_summary$cytoribo_id == cytoribo_id_curr,]
  
  # To use foldseek true hits
  foldseek_out_cytoribo_summary_curr_true_hit <- foldseek_out_cytoribo_curr[foldseek_out_cytoribo_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  foldseek_out_cytoribo_summary_curr <- data.frame(cytoribo_id_curr, human = 0, yeast = 0, trypanosoma = 0, giardia = 0)
  foldseek_out_cytoribo_summary_curr$human <- as.numeric("9606" %in% foldseek_out_cytoribo_curr$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$human <- foldseek_out_cytoribo_summary_curr$human + as.numeric("9606" %in% foldseek_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$yeast <- as.numeric("559292" %in% foldseek_out_cytoribo_curr$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$yeast <- foldseek_out_cytoribo_summary_curr$yeast + as.numeric("559292" %in% foldseek_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$trypanosoma <- as.numeric("185431" %in% foldseek_out_cytoribo_curr$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$trypanosoma <- foldseek_out_cytoribo_summary_curr$trypanosoma + as.numeric("185431" %in% foldseek_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$giardia <- as.numeric("5741" %in% foldseek_out_cytoribo_curr$target_taxid) * 0.5
  foldseek_out_cytoribo_summary_curr$giardia <- foldseek_out_cytoribo_summary_curr$giardia + as.numeric("5741" %in% foldseek_out_cytoribo_summary_curr_true_hit$target_taxid) * 0.5
  
  foldseek_cytoribo_species <- rbind(foldseek_cytoribo_species, foldseek_out_cytoribo_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(foldseek_cytoribo_species[,-1])
color_mat[color_mat == 0.5] <- 0

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- foldseek_cytoribo_species$cytoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(cytoribo_table_core$Subunit.Name))

# Plotting the tile heatmap with ggplot2
pdf("cytoribo_foldseek.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




## EggNOG cytoribo
# Idmapping for human and yeast cytoribo
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "cytoribo_human_uniprot_to_ensemblprotein_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "cytoribo_yeast_uniprot_to_kegg_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

eggnog_species_ids <- c("9606", "4932", "185431", "184922")

OG_cytoribo <- c()
for (i in 1:nrow(cytoribo_table_core)) {
  cytoribo_curr_row <- cytoribo_table_core[i,]
  
  # Map ids to eggnog
  tbr_geneid <- paste0("185431.", cytoribo_curr_row$Trypanosoma_brucei)
  giardia_geneid <- paste0("184922.", cytoribo_curr_row$Giardia_intestinalis)
  sce_geneid <- paste0("4932.", idmapping_yeast$To[idmapping_yeast$From %in% cytoribo_curr_row$Saccharomyces_cerevisiae])
  hsa_geneid <- paste0("9606.", idmapping_human$To[idmapping_human$From %in% cytoribo_curr_row$Homo_sapiens])
  
  cytoribo_curr_accessions <- c(hsa_geneid, sce_geneid, tbr_geneid, giardia_geneid)
  
  OG_curr_row <- root_OG_long %>% filter(protein_ids %in% cytoribo_curr_accessions)
  OG_curr_row$cytoribo_id <- cytoribo_curr_row$Subunit.Name
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "protein_ids", "cytoribo_id")]
  
  # Find missing entries
  missing_species <- eggnog_species_ids[!eggnog_species_ids %in% gsub("\\..*", "", OG_curr_row$protein_ids)]
  if (length(missing_species) > 0) {
    print(cytoribo_curr_row$Subunit.Name)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", protein_ids = missing_species, cytoribo_id = cytoribo_curr_row$Subunit.Name)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_cytoribo <- rbind(OG_cytoribo, OG_curr_row)
}

# Get taxids
OG_cytoribo$taxid <- gsub("\\..*", "", OG_cytoribo$protein_ids)

OG_cytoribo_summary <- OG_cytoribo %>% group_by(cytoribo_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_cytoribo_species <- c()
for (cytoribo_id_curr in unique(OG_cytoribo_summary$cytoribo_id)) {
  OG_cytoribo_summary_curr <- OG_cytoribo_summary[OG_cytoribo_summary$cytoribo_id == cytoribo_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 4932]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  giardia_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 184922]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  
  OG_cytoribo_species_curr <- data.frame(cytoribo_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Giardia = giardia_OG_id)
  OG_cytoribo_species <- rbind(OG_cytoribo_species, OG_cytoribo_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_cytoribo_species_bool <- t(apply(OG_cytoribo_species[,-1], 1, function(row) row == row[1]))
OG_cytoribo_species_bool[is.na(OG_cytoribo_species_bool)] <- FALSE
OG_cytoribo_species_colors <- ifelse(OG_cytoribo_species_bool, "#0047AB", "white")
OG_cytoribo_species_colors <- as.data.frame(OG_cytoribo_species_colors)

# Mark proteins missing from dataset as gray
OG_cytoribo_species_colors[which(OG_cytoribo_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_cytoribo_species_colors

color_df$row <- OG_cytoribo_species$cytoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(cytoribo_table_core$Subunit.Name))

# Plotting the tile heatmap with ggplot2
pdf("cytoribo_EggNOG_Eukaryota_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "EggNOG") +
  theme_minimal()
dev.off()



## OrthoDB cytoribo

# Idmapping for human and yeast cytoribo
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "cytoribo_human_uniprot_to_ensembl_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "cytoribo_yeast_uniprot_to_kegg_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

orthodb_species_ids <- c("185431_0", "559292_0", "9606_0", "5741_0")

OG_cytoribo <- c()
for (i in 1:nrow(cytoribo_table_core)) {
  cytoribo_curr_row <- cytoribo_table_core[i,]
  
  sce_geneid <- idmapping_yeast$To[idmapping_yeast$From %in% cytoribo_curr_row$Saccharomyces_cerevisiae]
  hsa_ensg_id <- idmapping_human$To[idmapping_human$From %in% cytoribo_curr_row$Homo_sapiens]
  
  tbr_gene_id <- idmapping_tbr$From[match(cytoribo_curr_row$Trypanosoma_brucei, idmapping_tbr$Entry)]
  
  # Convert accessions to RefSeq protein
  gla_refseq_protein_id <- gsub("^5741_", "", orthogroups_collapse_accessions_seprows_uniprot$accessions[match(cytoribo_curr_row$Giardia_intestinalis, orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession)])
  
  cytoribo_curr_accessions <- c(cytoribo_curr_row$Trypanosoma_brucei, cytoribo_curr_row$Giardia_intestinalis)
  cytoribo_curr_accessions <- cytoribo_curr_accessions[cytoribo_curr_accessions != ""]
  
  OG_curr_row <- odb_eukaryota_ogs_genes %>% filter(uniprot_id %in% cytoribo_curr_accessions | gene_id %in% hsa_ensg_id | synonym %in% sce_geneid | synonym %in% tbr_gene_id | seq_id %in% gla_refseq_protein_id)
  OG_curr_row$cytoribo_id <- cytoribo_curr_row$Subunit.Name
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "species_id", "cytoribo_id")]
  
  # Find missing entries
  missing_species <- orthodb_species_ids[!orthodb_species_ids %in% OG_curr_row$species_id]
  if (length(missing_species) > 0) {
    print(cytoribo_curr_row$Subunit.Name)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", species_id = missing_species, cytoribo_id = cytoribo_curr_row$Subunit.Name)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_cytoribo <- rbind(OG_cytoribo, OG_curr_row)
}

# Get taxids
OG_cytoribo$taxid <- gsub("_.*", "", OG_cytoribo$species_id)

OG_cytoribo_summary <- OG_cytoribo %>% group_by(cytoribo_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_cytoribo_species <- c()
for (cytoribo_id_curr in unique(OG_cytoribo_summary$cytoribo_id)) {
  OG_cytoribo_summary_curr <- OG_cytoribo_summary[OG_cytoribo_summary$cytoribo_id == cytoribo_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  giardia_OG_id <- OG_cytoribo_summary_curr$OG_id[OG_cytoribo_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  
  OG_cytoribo_species_curr <- data.frame(cytoribo_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Giardia = giardia_OG_id)
  OG_cytoribo_species <- rbind(OG_cytoribo_species, OG_cytoribo_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_cytoribo_species_bool <- t(apply(OG_cytoribo_species[,-1], 1, function(row) row == row[1]))
OG_cytoribo_species_bool[is.na(OG_cytoribo_species_bool)] <- FALSE
OG_cytoribo_species_colors <- ifelse(OG_cytoribo_species_bool, "#0047AB", "white")
OG_cytoribo_species_colors <- as.data.frame(OG_cytoribo_species_colors)

# Mark proteins missing from dataset as gray
OG_cytoribo_species_colors[which(OG_cytoribo_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_cytoribo_species_colors

color_df$row <- OG_cytoribo_species$cytoribo_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(cytoribo_table_core$Subunit.Name))


# Plotting the tile heatmap with ggplot2
pdf("cytoribo_OrthoDB_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "OrthoDB") +
  theme_minimal()
dev.off()



# Read in ISC gene sets
ISC_table_core <- read.csv(here("data", "annotation", "pathways", "ISC_pathway_annotation.csv"))
ISC_table_core <- ISC_table_core %>% rowwise() %>% filter(!Human_Symbol %in% c("BOLA1", "BOLA3"))


# OG data ISC
OG_ISC <- c()
for (i in 1:nrow(ISC_table_core)) {
  ISC_curr_row <- ISC_table_core[i,]
  
  ISC_curr_accessions <- c(ISC_curr_row$Homo_sapiens, ISC_curr_row$Saccharomyces_cerevisiae, ISC_curr_row$Trypanosoma_brucei, ISC_curr_row$Arabidopsis_thaliana_mito, ISC_curr_row$Arabidopsis_thaliana_plastid)

  ISC_curr_accessions <- unlist(strsplit(ISC_curr_accessions, ","))
  
  # Only include homologs from structure
  OG_curr_row <- orthogroups_collapse_accessions_seprows_uniprot[which(orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession %in% ISC_curr_accessions),]
  OG_curr_row$ISC_id <- ISC_curr_row$Human_Symbol
  
  OG_ISC <- rbind(OG_ISC, OG_curr_row)
}

# Get taxids
OG_ISC$taxid <- gsub("_.*", "", OG_ISC$accessions)

phrogs_long_select_human <- phrogs_long %>% filter(protein_id %in% OG_ISC$accessions) %>% filter(taxid == "9606")
phrogs_long_select_ISC <- phrogs_long %>% filter(PROG_id %in% phrogs_long_select_human$PROG_id)

if (BOOL_USE_PHROG) {
  OG_ISC <- OG_ISC %>% filter(accessions %in% phrogs_long_select_ISC$protein_id)
} else {
  OG_ISC <- OG_ISC %>% group_by(ISC_id) %>% filter(OG_id %in% OG_id[taxid == "9606"]) %>% ungroup()
}

if (BOOL_USE_ABSENSE) {
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows, phrogs_long_select_human, by.x = c("human_accessions", "OG_id"), by.y = c("protein_id", "OG_id"))
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows_human, OG_ISC, by.x = c("human_accessions", "OG_id"), by.y = c("accessions", "OG_id"))
  homology_power_per_species_curr <- homology_power_per_species %>% filter(PhROG_id %in% homology_power_agg_seprows_human$PhROG_id) %>% filter(taxid %in% c("559292", "3702", "185431"))
  homology_power_per_species_curr$ISC_id <- homology_power_agg_seprows_human$ISC_id[match(homology_power_per_species_curr$PhROG_id, homology_power_agg_seprows_human$PhROG_id)]
  homology_power_per_species_human_underpowered <- homology_power_per_species_curr %>% filter(probability_of_detection < 0.95)
  
  if (nrow(homology_power_per_species_human_underpowered) > 0) {
    OG_ISC_underpowered <- data.frame(OG_id = "Underpowered", accessions = NA, uniprot_accession = NA, ISC_id = homology_power_per_species_human_underpowered$ISC_id, taxid = homology_power_per_species_human_underpowered$taxid)
    OG_ISC <- rbind(OG_ISC, OG_ISC_underpowered)
    OG_ISC <- OG_ISC %>% group_by(ISC_id) %>% filter(!duplicated(taxid)) %>% ungroup()
  }
}

# Mark Ath plastid
ath_plastid_ids <- unlist(strsplit(ISC_table_core$Arabidopsis_thaliana_plastid, ","))
OG_ISC$taxid[OG_ISC$uniprot_accession %in% ath_plastid_ids] <- paste0(OG_ISC$taxid[OG_ISC$uniprot_accession %in% ath_plastid_ids], "_plastid")

OG_ISC_summary <- OG_ISC %>% group_by(ISC_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_ISC_species <- c()
for (ISC_id_curr in unique(OG_ISC_summary$ISC_id)) {
  OG_ISC_summary_curr <- OG_ISC_summary[OG_ISC_summary$ISC_id == ISC_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  arabidopsis_plastid_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == "3702_plastid"]
  if (length(arabidopsis_plastid_OG_id) == 0) {
    arabidopsis_plastid_OG_id <- NA
  }
  
  OG_ISC_species_curr <- data.frame(ISC_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Arabidopsis = arabidopsis_OG_id, Arabidopsis_plastid = arabidopsis_plastid_OG_id)
  OG_ISC_species <- rbind(OG_ISC_species, OG_ISC_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Assign colors
if (BOOL_USE_ABSENSE) {
  OG_ISC_species_colors <- t(apply(OG_ISC_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", row == "Underpowered" ~ "gray", TRUE ~ "white")))
} else {
  OG_ISC_species_colors <- t(apply(OG_ISC_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", TRUE ~ "white")))
}
colnames(OG_ISC_species_colors) <- colnames(OG_ISC_species[,-1])
OG_ISC_species_colors <- as.data.frame(OG_ISC_species_colors)

# Reshape the color data into a dataframe for plotting
color_df <- OG_ISC_species_colors

color_df$row <- OG_ISC_species$ISC_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(ISC_table_core$Human_Symbol))

if (BOOL_USE_PHROG) {
  dataset_description <- "PhROG"
} else {
  dataset_description <- "OG"
}

# Plotting the tile heatmap with ggplot2
pdf(paste0("ISC_", dataset_description, "_heatmap.pdf"), width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = dataset_description) +
  theme_minimal()
dev.off()




# Blast analysis
blast_out <- blast_out_BBH
title_string <- "BLAST BBH expect<1 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(blast_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
blast_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
blast_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(blast_out$sseqid, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

# Separate human rows to capture paralogs
ISC_table_core_seprows_human <- ISC_table_core %>% separate_rows(Homo_sapiens, sep=",")

blast_out_ISC <- blast_out[which(blast_out$query_uniprot_accession %in% ISC_table_core_seprows_human$Homo_sapiens),]
blast_out_ISC$ISC_id <- ISC_table_core_seprows_human$Human_Symbol[match(blast_out_ISC$query_uniprot_accession, ISC_table_core_seprows_human$Homo_sapiens)]
blast_out_ISC <- blast_out_ISC %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% unlist(strsplit(unlist(ISC_table_core[which(ISC_table_core$Human_Symbol == ISC_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Arabidopsis_thaliana_mito", "Arabidopsis_thaliana_plastid")]), ",")))
blast_out_ISC_summary <- blast_out_ISC %>% group_by(ISC_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Split into 1 column per species
blast_ISC_species <- c()
for (ISC_id_curr in unique(blast_out_ISC_summary$ISC_id)) {
  blast_out_ISC_curr <- blast_out_ISC_summary[blast_out_ISC_summary$ISC_id == ISC_id_curr,]
  
  # To use blast true hits
  blast_out_ISC_summary_curr_true_hit <- blast_out_ISC_curr[blast_out_ISC_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  blast_out_ISC_summary_curr <- data.frame(ISC_id_curr, human = 0, yeast = 0, trypanosoma = 0, arabidopsis = 0, arabidopsis_plastid = 0)
  blast_out_ISC_summary_curr$human <- as.numeric("9606" %in% blast_out_ISC_curr$target_taxid) * 0.5
  blast_out_ISC_summary_curr$human <- blast_out_ISC_summary_curr$human + as.numeric("9606" %in% blast_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_ISC_summary_curr$yeast <- as.numeric("559292" %in% blast_out_ISC_curr$target_taxid) * 0.5
  blast_out_ISC_summary_curr$yeast <- blast_out_ISC_summary_curr$yeast + as.numeric("559292" %in% blast_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_ISC_summary_curr$trypanosoma <- as.numeric("185431" %in% blast_out_ISC_curr$target_taxid) * 0.5
  blast_out_ISC_summary_curr$trypanosoma <- blast_out_ISC_summary_curr$trypanosoma + as.numeric("185431" %in% blast_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_ISC_summary_curr$arabidopsis <- as.numeric("3702" %in% blast_out_ISC_curr$target_taxid) * 0.5
  blast_out_ISC_summary_curr$arabidopsis <- blast_out_ISC_summary_curr$arabidopsis + as.numeric("3702" %in% blast_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_ISC_summary_curr$arabidopsis_plastid <- as.numeric("3702_plastid" %in% blast_out_ISC_curr$target_taxid) * 0.5
  blast_out_ISC_summary_curr$arabidopsis_plastid <- blast_out_ISC_summary_curr$arabidopsis_plastid + as.numeric("3702_plastid" %in% blast_out_ISC_summary_curr_true_hit$target_taxid) * 0.5

  blast_ISC_species <- rbind(blast_ISC_species, blast_out_ISC_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(blast_ISC_species[,-1])
color_mat[color_mat == 0.5] <- 0

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- blast_ISC_species$ISC_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(ISC_table_core$Human_Symbol))

# Plotting the tile heatmap with ggplot2
pdf("ISC_blast.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "BLAST") +
  theme_minimal()
# theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




# plmsearch analysis
plmsearch_out <- plmsearch_out_BBH
title_string <- "plmsearch BBH similarity score > 0.3 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(plmsearch_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
plmsearch_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
plmsearch_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(plmsearch_out$target, orthogroups_collapse_accessions_seprows_uniprot$accessions)]

# Separate human rows to capture paralogs
ISC_table_core_seprows_human <- ISC_table_core %>% separate_rows(Homo_sapiens, sep=",")

plmsearch_out_ISC <- plmsearch_out[which(plmsearch_out$query_uniprot_accession %in% ISC_table_core_seprows_human$Homo_sapiens),]
plmsearch_out_ISC$ISC_id <- ISC_table_core_seprows_human$Human_Symbol[match(plmsearch_out_ISC$query_uniprot_accession, ISC_table_core_seprows_human$Homo_sapiens)]
plmsearch_out_ISC <- plmsearch_out_ISC %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% unlist(strsplit(unlist(ISC_table_core[which(ISC_table_core$Human_Symbol == ISC_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Arabidopsis_thaliana_mito", "Arabidopsis_thaliana_plastid")]), ",")))
plmsearch_out_ISC_summary <- plmsearch_out_ISC %>% group_by(ISC_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Split into 1 column per species
plmsearch_ISC_species <- c()
for (ISC_id_curr in unique(plmsearch_out_ISC_summary$ISC_id)) {
  plmsearch_out_ISC_curr <- plmsearch_out_ISC_summary[plmsearch_out_ISC_summary$ISC_id == ISC_id_curr,]
  
  # To use plmsearch true hits
  plmsearch_out_ISC_summary_curr_true_hit <- plmsearch_out_ISC_curr[plmsearch_out_ISC_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  plmsearch_out_ISC_summary_curr <- data.frame(ISC_id_curr, human = 0, yeast = 0, trypanosoma = 0, arabidopsis = 0, arabidopsis_plastid = 0)
  plmsearch_out_ISC_summary_curr$human <- as.numeric("9606" %in% plmsearch_out_ISC_curr$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$human <- plmsearch_out_ISC_summary_curr$human + as.numeric("9606" %in% plmsearch_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$yeast <- as.numeric("559292" %in% plmsearch_out_ISC_curr$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$yeast <- plmsearch_out_ISC_summary_curr$yeast + as.numeric("559292" %in% plmsearch_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$trypanosoma <- as.numeric("185431" %in% plmsearch_out_ISC_curr$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$trypanosoma <- plmsearch_out_ISC_summary_curr$trypanosoma + as.numeric("185431" %in% plmsearch_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$arabidopsis <- as.numeric("3702" %in% plmsearch_out_ISC_curr$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$arabidopsis <- plmsearch_out_ISC_summary_curr$arabidopsis + as.numeric("3702" %in% plmsearch_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$arabidopsis_plastid <- as.numeric("3702_plastid" %in% plmsearch_out_ISC_curr$target_taxid) * 0.5
  plmsearch_out_ISC_summary_curr$arabidopsis_plastid <- plmsearch_out_ISC_summary_curr$arabidopsis_plastid + as.numeric("3702_plastid" %in% plmsearch_out_ISC_summary_curr_true_hit$target_taxid) * 0.5

  plmsearch_ISC_species <- rbind(plmsearch_ISC_species, plmsearch_out_ISC_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(plmsearch_ISC_species[,-1])
color_mat[color_mat == 0.5] <- 0

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- plmsearch_ISC_species$ISC_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(ISC_table_core$Human_Symbol))

# Plotting the tile heatmap with ggplot2
pdf("ISC_plmsearch.plmalign.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "plmsearch") +
  theme_minimal()
dev.off()




# foldseek analysis
foldseek_out <- foldseek_out_BBH
title_string <- "foldseek BBH expect<1 with human homolog"

# Separate human rows to capture paralogs
ISC_table_core_seprows_human <- ISC_table_core %>% separate_rows(Homo_sapiens, sep=",")

foldseek_out_ISC <- foldseek_out[which(foldseek_out$query %in% ISC_table_core_seprows_human$Homo_sapiens),]
foldseek_out_ISC$ISC_id <- ISC_table_core_seprows_human$Human_Symbol[match(foldseek_out_ISC$query, ISC_table_core_seprows_human$Homo_sapiens)]
foldseek_out_ISC <- foldseek_out_ISC %>% rowwise() %>% mutate(true_top_hit = target %in% c(unlist(strsplit(unlist(ISC_table_core[which(ISC_table_core$Human_Symbol == ISC_id), c("Homo_sapiens", "Saccharomyces_cerevisiae", "Trypanosoma_brucei", "Arabidopsis_thaliana_mito", "Arabidopsis_thaliana_plastid")]), ",")), ath_idmapping$Entry[which(ath_idmapping$From %in% ath_idmapping$From[match(unlist(strsplit(ISC_table_core$Arabidopsis_thaliana_mito[which(ISC_table_core$Human_Symbol == ISC_id)], ",")), ath_idmapping$Entry)])], ath_idmapping$Entry[which(ath_idmapping$From %in% ath_idmapping$From[match(unlist(strsplit(ISC_table_core$Arabidopsis_thaliana_plastid[which(ISC_table_core$Human_Symbol == ISC_id)], ",")), ath_idmapping$Entry)])]))
foldseek_out_ISC_summary <- foldseek_out_ISC %>% group_by(ISC_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Split into 1 column per species
foldseek_ISC_species <- c()
for (ISC_id_curr in unique(foldseek_out_ISC_summary$ISC_id)) {
  foldseek_out_ISC_curr <- foldseek_out_ISC_summary[foldseek_out_ISC_summary$ISC_id == ISC_id_curr,]
  
  # To use foldseek true hits
  foldseek_out_ISC_summary_curr_true_hit <- foldseek_out_ISC_curr[foldseek_out_ISC_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  foldseek_out_ISC_summary_curr <- data.frame(ISC_id_curr, human = 0, yeast = 0, trypanosoma = 0, arabidopsis = 0, arabidopsis_plastid = 0)
  foldseek_out_ISC_summary_curr$human <- as.numeric("9606" %in% foldseek_out_ISC_curr$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$human <- foldseek_out_ISC_summary_curr$human + as.numeric("9606" %in% foldseek_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$yeast <- as.numeric("559292" %in% foldseek_out_ISC_curr$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$yeast <- foldseek_out_ISC_summary_curr$yeast + as.numeric("559292" %in% foldseek_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$trypanosoma <- as.numeric("185431" %in% foldseek_out_ISC_curr$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$trypanosoma <- foldseek_out_ISC_summary_curr$trypanosoma + as.numeric("185431" %in% foldseek_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$arabidopsis <- as.numeric("3702" %in% foldseek_out_ISC_curr$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$arabidopsis <- foldseek_out_ISC_summary_curr$arabidopsis + as.numeric("3702" %in% foldseek_out_ISC_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$arabidopsis_plastid <- as.numeric("3702_plastid" %in% foldseek_out_ISC_curr$target_taxid) * 0.5
  foldseek_out_ISC_summary_curr$arabidopsis_plastid <- foldseek_out_ISC_summary_curr$arabidopsis_plastid + as.numeric("3702_plastid" %in% foldseek_out_ISC_summary_curr_true_hit$target_taxid) * 0.5

  foldseek_ISC_species <- rbind(foldseek_ISC_species, foldseek_out_ISC_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(foldseek_ISC_species[,-1])
color_mat[color_mat == 0.5] <- 0

# Binary color
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- foldseek_ISC_species$ISC_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(ISC_table_core$Human_Symbol))

# Plotting the tile heatmap with ggplot2
pdf("ISC_foldseek.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "foldseek") +
  theme_minimal()
dev.off()




## EggNOG ISC
# Idmapping for human and yeast ISC
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "ISC_human_uniprot_to_ensemblprotein_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "ISC_yeast_uniprot_to_kegg_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

eggnog_species_ids <- c("9606", "4932", "3702", "185431")

OG_ISC <- c()
for (i in 1:nrow(ISC_table_core)) {
  ISC_curr_row <- ISC_table_core[i,]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_ids <- unlist(strsplit(ISC_curr_row$Arabidopsis_thaliana_mito, ","))
  ath_geneid_curr <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- paste0("3702.", ath_idmapping$Entry[which(ath_idmapping$From %in% ath_geneid_curr)])
  
  # Map ids to eggnog
  ath_geneid <- paste0("3702.", unlist(strsplit(ath_ids, ",")))
  tbr_geneid <- paste0("185431.", unlist(strsplit(ISC_curr_row$Trypanosoma_brucei, ",")))
  sce_geneid <- paste0("4932.", idmapping_yeast$To[idmapping_yeast$From %in% unlist(strsplit(ISC_curr_row$Saccharomyces_cerevisiae, ","))])
  hsa_geneid <- paste0("9606.", idmapping_human$To[idmapping_human$From %in% unlist(strsplit(ISC_curr_row$Homo_sapiens, ","))])
  
  ISC_curr_accessions <- c(hsa_geneid, sce_geneid, ath_geneid, tbr_geneid, ath_uniprot_ids)
  
  OG_curr_row <- root_OG_long %>% filter(protein_ids %in% ISC_curr_accessions)
  
  ## Get plastid proteins
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_ids <- unlist(strsplit(ISC_curr_row$Arabidopsis_thaliana_plastid, ","))
  ath_geneid_curr <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- paste0("3702.", ath_idmapping$Entry[which(ath_idmapping$From %in% ath_geneid_curr)])
  ath_geneid <- paste0("3702.", unlist(strsplit(ath_ids, ",")))
  OG_curr_row_plastid <- root_OG_long %>% filter(protein_ids %in% c(ath_geneid, ath_uniprot_ids))
  OG_curr_row_plastid$protein_ids <- gsub("3702.", "3702_plastid.", OG_curr_row_plastid$protein_ids, fixed=TRUE)
  OG_curr_row <- rbind(OG_curr_row, OG_curr_row_plastid)
  
  OG_curr_row$ISC_id <- ISC_curr_row$Human_Symbol
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "protein_ids", "ISC_id")]
  
  # Find missing entries
  missing_species <- eggnog_species_ids[!eggnog_species_ids %in% gsub("\\..*", "", OG_curr_row$protein_ids)]
  if (length(missing_species) > 0) {
    print(ISC_curr_row$Human_Symbol)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", protein_ids = missing_species, ISC_id = ISC_curr_row$Human_Symbol)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_ISC <- rbind(OG_ISC, OG_curr_row)
}

# Get taxids
OG_ISC$taxid <- gsub("\\..*", "", OG_ISC$protein_ids)

OG_ISC_summary <- OG_ISC %>% group_by(ISC_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_ISC_species <- c()
for (ISC_id_curr in unique(OG_ISC_summary$ISC_id)) {
  OG_ISC_summary_curr <- OG_ISC_summary[OG_ISC_summary$ISC_id == ISC_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 4932]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  arabidopsis_plastid_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == "3702_plastid"]
  if (length(arabidopsis_plastid_OG_id) == 0) {
    arabidopsis_plastid_OG_id <- NA
  }
  
  OG_ISC_species_curr <- data.frame(ISC_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Arabidopsis = arabidopsis_OG_id, Arabidopsis_plastid = arabidopsis_plastid_OG_id)
  OG_ISC_species <- rbind(OG_ISC_species, OG_ISC_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_ISC_species_bool <- t(apply(OG_ISC_species[,-1], 1, function(row) row == row[1]))
OG_ISC_species_bool[is.na(OG_ISC_species_bool)] <- FALSE
OG_ISC_species_colors <- ifelse(OG_ISC_species_bool, "#0047AB", "white")
OG_ISC_species_colors <- as.data.frame(OG_ISC_species_colors)

# Mark proteins missing from dataset as gray
OG_ISC_species_colors[which(OG_ISC_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_ISC_species_colors

color_df$row <- OG_ISC_species$ISC_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(ISC_table_core$Human_Symbol))

# Plotting the tile heatmap with ggplot2
pdf("ISC_EggNOG_Eukaryota_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "EggNOG") +
  theme_minimal()
dev.off()



# OrthoDB ISC
# Idmapping for human and yeast ISC
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "ISC_human_uniprot_to_ensembl_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

idmapping_yeast <- read.table(here("data", "annotation", "pathways", "idmapping", "ISC_yeast_uniprot_to_kegg_idmapping_2025_07_08.tsv"), sep="\t", header=TRUE)
idmapping_yeast$To <- gsub("sce:", "", idmapping_yeast$To, fixed=TRUE)

orthodb_species_ids <- c("185431_0", "3702_0", "559292_0", "9606_0")

OG_ISC <- c()
for (i in 1:nrow(ISC_table_core)) {
  ISC_curr_row <- ISC_table_core[i,]
  
  sce_geneid <- idmapping_yeast$To[idmapping_yeast$From %in% unlist(strsplit(ISC_curr_row$Saccharomyces_cerevisiae, ","))]
  hsa_ensg_id <- idmapping_human$To[idmapping_human$From %in% unlist(strsplit(ISC_curr_row$Homo_sapiens, ","))]
  
  tbr_gene_id <- idmapping_tbr$From[match(unlist(strsplit(ISC_curr_row$Trypanosoma_brucei, ",")), idmapping_tbr$Entry)]
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_ids <- unlist(strsplit(ISC_curr_row$Arabidopsis_thaliana_mito, ","))
  ath_geneid <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From %in% ath_geneid)]
  idmapping_ath_curr <- idmapping_ath %>% filter(KEGG %in% ath_geneid | TAIR %in% ath_geneid | Araport %in% ath_geneid)
  ath_uniprot_ids <- c(ath_uniprot_ids, idmapping_ath_curr$From)
  
  ISC_curr_accessions <- c(ISC_curr_row$Trypanosoma_brucei, ath_ids, ath_uniprot_ids)
  ISC_curr_accessions <- ISC_curr_accessions[ISC_curr_accessions != ""]
  
  OG_curr_row <- odb_eukaryota_ogs_genes %>% filter(uniprot_id %in% ISC_curr_accessions | gene_id %in% hsa_ensg_id | synonym %in% sce_geneid | gene_id %in% ath_geneid | synonym %in% tbr_gene_id)
  
  ## Get Arabidopsis plastid proteins
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_ids <- unlist(strsplit(ISC_curr_row$Arabidopsis_thaliana_plastid, ","))
  ath_geneid <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From %in% ath_geneid)]
  idmapping_ath_curr <- idmapping_ath %>% filter(KEGG %in% ath_geneid | TAIR %in% ath_geneid | Araport %in% ath_geneid)
  ath_uniprot_ids <- c(ath_uniprot_ids, idmapping_ath_curr$From)
  
  ISC_curr_accessions <- c(ath_ids, ath_uniprot_ids)
  ISC_curr_accessions <- ISC_curr_accessions[ISC_curr_accessions != ""]
  OG_curr_row_plastid <- odb_eukaryota_ogs_genes %>% filter(uniprot_id %in% ISC_curr_accessions | gene_id %in% ath_geneid)
  OG_curr_row_plastid$species_id <- gsub("3702_0", "3702plastid_0", OG_curr_row_plastid$species_id, fixed=TRUE)
  OG_curr_row <- rbind(OG_curr_row, OG_curr_row_plastid)
  
  OG_curr_row$ISC_id <- ISC_curr_row$Human_Symbol
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "species_id", "ISC_id")]
  
  # Find missing entries
  missing_species <- orthodb_species_ids[!orthodb_species_ids %in% OG_curr_row$species_id]
  if (length(missing_species) > 0) {
    print(ISC_curr_row$Human_Symbol)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", species_id = missing_species, ISC_id = ISC_curr_row$Human_Symbol)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_ISC <- rbind(OG_ISC, OG_curr_row)
}

# Get taxids
OG_ISC$taxid <- gsub("_.*", "", OG_ISC$species_id)

OG_ISC_summary <- OG_ISC %>% group_by(ISC_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_ISC_species <- c()
for (ISC_id_curr in unique(OG_ISC_summary$ISC_id)) {
  OG_ISC_summary_curr <- OG_ISC_summary[OG_ISC_summary$ISC_id == ISC_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  trypanosoma_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 185431]
  if (length(trypanosoma_OG_id) == 0) {
    trypanosoma_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  arabidopsis_plastid_OG_id <- OG_ISC_summary_curr$OG_id[OG_ISC_summary_curr$taxid == "3702plastid"]
  if (length(arabidopsis_plastid_OG_id) == 0) {
    arabidopsis_plastid_OG_id <- NA
  }
  
  OG_ISC_species_curr <- data.frame(ISC_id_curr, Human = human_OG_id, Yeast = yeast_OG_id, Trypanosoma = trypanosoma_OG_id, Arabidopsis = arabidopsis_OG_id, Arabidopsis_plastid = arabidopsis_plastid_OG_id)
  OG_ISC_species <- rbind(OG_ISC_species, OG_ISC_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_ISC_species_bool <- t(apply(OG_ISC_species[,-1], 1, function(row) row == row[1]))
OG_ISC_species_bool[is.na(OG_ISC_species_bool)] <- FALSE
OG_ISC_species_colors <- ifelse(OG_ISC_species_bool, "#0047AB", "white")
OG_ISC_species_colors <- as.data.frame(OG_ISC_species_colors)

# Mark proteins missing from dataset as gray
OG_ISC_species_colors[which(OG_ISC_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_ISC_species_colors

color_df$row <- OG_ISC_species$ISC_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(ISC_table_core$Human_Symbol))


# Plotting the tile heatmap with ggplot2
pdf("ISC_OrthoDB_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "OrthoDB") +
  theme_minimal()
dev.off()



# Read in Complex I data
complexI_table <- read.csv(here("data", "annotation", "pathways", "complexI_structure_annotation.csv"))
complexI_table_core <- complexI_table %>% rowwise() %>% filter(!Homo_sapiens %in% c("", "missing_from_uniprot") & !Tetrahymena_thermophila %in% c("", "missing_from_uniprot") & !Arabidopsis_thaliana %in% c("", "missing_from_uniprot"))
complex_table_core <- complexI_table_core

# OG data complex I
OG_complex <- c()
for (i in 1:nrow(complex_table_core)) {
  complex_curr_row <- complex_table_core[i,]
  
  # Get all Arabidopsis paralogs
  ath_ids <- unlist(strsplit(complex_curr_row$Arabidopsis_thaliana, ","))
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid)]
  
  complex_curr_accessions <- c(complex_curr_row$Homo_sapiens, complex_curr_row$Tetrahymena_thermophila, ath_ids, complex_curr_row$Spinach_chloroplast, ath_uniprot_ids)
  
  # Only include homologs from structure
  OG_curr_row <- orthogroups_collapse_accessions_seprows_uniprot[which(orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession %in% complex_curr_accessions),]
  
  # Add tetrahymena mtDNA proteins
  tetrahymena_accession <- paste0("312017_", complex_curr_row$Tetrahymena_thermophila)
  if (tetrahymena_accession %in% phrogs_long$protein_id) {
    OG_curr_row_tetrahymena <- data.frame(OG_id = phrogs_long$OG_id[phrogs_long$protein_id == tetrahymena_accession], accessions = tetrahymena_accession, uniprot_accession = tetrahymena_accession)
    OG_curr_row <- rbind(OG_curr_row, OG_curr_row_tetrahymena)
  }
  
  OG_curr_row$complex_id <- complex_curr_row$Name
  
  # Include homologs from other species
  OG_curr_all_members <- orthogroups_collapse_accessions_seprows_uniprot[orthogroups_collapse_accessions_seprows_uniprot$OG_id %in% OG_curr_row$OG_id,]
  # OG_curr_all_members <- OG_curr_all_members[gsub("_.*", "", OG_curr_all_members$accessions) %in% c("9606", "559292", "5702", "3702"),]
  yeast_index <- which(gsub("_.*", "", OG_curr_all_members$accessions) %in% c("559292"))
  if (length(yeast_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[yeast_index,]
    OG_curr_all_members$complex_id <- complex_curr_row$Name
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_complex <- rbind(OG_complex, OG_curr_row)
}

# Get taxids
OG_complex$taxid <- gsub("_.*", "", OG_complex$accessions)

phrogs_long_select_human <- phrogs_long %>% filter(protein_id %in% OG_complex$accessions) %>% filter(taxid == "9606")
phrogs_long_select_complex <- phrogs_long %>% filter(PROG_id %in% phrogs_long_select_human$PROG_id)

if (BOOL_USE_PHROG) {
  OG_complex <- OG_complex %>% filter(accessions %in% phrogs_long_select_complex$protein_id)
} else {
  OG_complex <- OG_complex %>% group_by(complex_id) %>% filter(OG_id %in% OG_id[taxid == "9606"]) %>% ungroup()
}

if (BOOL_USE_ABSENSE) {
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows, phrogs_long_select_human, by.x = c("human_accessions", "OG_id"), by.y = c("protein_id", "OG_id"))
  homology_power_agg_seprows_human <- merge(homology_power_agg_seprows_human, OG_complex, by.x = c("human_accessions", "OG_id"), by.y = c("accessions", "OG_id"))
  homology_power_per_species_curr <- homology_power_per_species %>% filter(PhROG_id %in% homology_power_agg_seprows_human$PhROG_id) %>% filter(taxid %in% c("559292", "3702", "312017", "5741", "3562"))
  homology_power_per_species_curr$complex_id <- homology_power_agg_seprows_human$complex_id[match(homology_power_per_species_curr$PhROG_id, homology_power_agg_seprows_human$PhROG_id)]
  homology_power_per_species_human_underpowered <- homology_power_per_species_curr %>% filter(probability_of_detection < 0.95)
  OG_complex_underpowered <- data.frame(OG_id = "Underpowered", accessions = NA, uniprot_accession = NA, complex_id = homology_power_per_species_human_underpowered$complex_id, taxid = homology_power_per_species_human_underpowered$taxid)
  OG_complex <- rbind(OG_complex, OG_complex_underpowered)
  OG_complex <- OG_complex %>% group_by(complex_id) %>% filter(!duplicated(taxid)) %>% ungroup()
}

OG_complex_summary <- OG_complex %>% group_by(complex_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_complex_species <- c()

for (complex_id_curr in unique(OG_complex_summary$complex_id)) {
  OG_complex_summary_curr <- OG_complex_summary[OG_complex_summary$complex_id == complex_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  tetrahymena_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 312017]
  if (length(tetrahymena_OG_id) == 0) {
    tetrahymena_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  giardia_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  spinach_chloroplast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_OG_id) == 0) {
    spinach_chloroplast_OG_id <- NA
  }
  
  OG_complex_species_curr <- data.frame(complex_id_curr, `9606` = human_OG_id, `559292` = yeast_OG_id, `312017` = tetrahymena_OG_id, `3702` = arabidopsis_OG_id, `5741` = giardia_OG_id, `3562` = spinach_chloroplast_OG_id, check.names = FALSE)
  OG_complex_species <- rbind(OG_complex_species, OG_complex_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Assign colors
if (BOOL_USE_ABSENSE) {
  OG_complex_species_colors <- t(apply(OG_complex_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", row == "Underpowered" ~ "gray", TRUE ~ "white")))
} else {
  OG_complex_species_colors <- t(apply(OG_complex_species[,-1], 1, function(row) case_when(row == row[1] ~ "#0047AB", TRUE ~ "white")))
}
colnames(OG_complex_species_colors) <- colnames(OG_complex_species[,-1])
OG_complex_species_colors <- as.data.frame(OG_complex_species_colors)

# Reshape the color data into a dataframe for plotting
color_df <- OG_complex_species_colors

color_df$row <- OG_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$Name))

if (BOOL_USE_PHROG) {
  dataset_description <- "PhROG"
} else {
  dataset_description <- "OG"
}

# Plotting the tile heatmap with ggplot2
pdf(paste0("complexI_", dataset_description, "_heatmap.pdf"), width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = dataset_description) +
  theme_minimal()
dev.off()




# Blast analysis
blast_out <- blast_out_BBH
title_string <- "BLAST BBH expect<1 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(blast_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
blast_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
blast_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(blast_out$sseqid, orthogroups_collapse_accessions_seprows_uniprot$accessions)]
# Keep mtDNA accessions (starting with AAK) for Tetrahymena
blast_out$target_uniprot_accession[blast_out$target_taxid == 312017 & grepl("312017_AAK", blast_out$sseqid, fixed=TRUE)] <- blast_out$target_uniprot_entry[blast_out$target_taxid == 312017 & grepl("312017_AAK", blast_out$sseqid, fixed=TRUE)]

blast_out_complex <- blast_out[which(blast_out$query_uniprot_accession %in% complex_table_core$Homo_sapiens),]
blast_out_complex$complex_id <- complex_table_core$Name[match(blast_out_complex$query_uniprot_accession, complex_table_core$Homo_sapiens)]
blast_out_complex <- blast_out_complex %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% unlist(strsplit(unlist(complex_table_core[which(complex_table_core$Name == complex_id), c("Homo_sapiens", "Tetrahymena_thermophila", "Arabidopsis_thaliana", "Spinach_chloroplast")]), ",")))
blast_out_complex_summary <- blast_out_complex %>% group_by(complex_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Yeast and Giardia hits as true since has been lost
blast_out_complex_summary$any_true_top_hit[blast_out_complex_summary$target_taxid == "559292"] <- TRUE 
blast_out_complex_summary$any_true_top_hit[blast_out_complex_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
blast_complex_species <- c()
for (complex_id_curr in unique(blast_out_complex_summary$complex_id)) {
  blast_out_complex_curr <- blast_out_complex_summary[blast_out_complex_summary$complex_id == complex_id_curr,]
  
  # To use blast true hits
  blast_out_complex_summary_curr_true_hit <- blast_out_complex_curr[blast_out_complex_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  blast_out_complex_summary_curr <- data.frame(complex_id_curr, human = 0, yeast = 0, tetrahymena = 0, arabidopsis = 0, giardia = 0, spinach_chloroplast = 0)
  
  # Binary
  blast_out_complex_summary_curr$human <- as.numeric("9606" %in% blast_out_complex_summary_curr_true_hit$target_taxid)
  blast_out_complex_summary_curr$yeast <- as.numeric("559292" %in% blast_out_complex_summary_curr_true_hit$target_taxid)
  blast_out_complex_summary_curr$tetrahymena <- as.numeric("312017" %in% blast_out_complex_summary_curr_true_hit$target_taxid)
  blast_out_complex_summary_curr$arabidopsis <- as.numeric("3702" %in% blast_out_complex_summary_curr_true_hit$target_taxid)
  blast_out_complex_summary_curr$giardia <- as.numeric("5741" %in% blast_out_complex_summary_curr_true_hit$target_taxid)
  blast_out_complex_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% blast_out_complex_summary_curr_true_hit$target_taxid)
  
  # Tri color
  blast_out_complex_summary_curr$human <- as.numeric("9606" %in% blast_out_complex_curr$target_taxid) * 0.5
  blast_out_complex_summary_curr$human <- blast_out_complex_summary_curr$human + as.numeric("9606" %in% blast_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_complex_summary_curr$tetrahymena <- as.numeric("312017" %in% blast_out_complex_curr$target_taxid) * 0.5
  blast_out_complex_summary_curr$tetrahymena <- blast_out_complex_summary_curr$tetrahymena + as.numeric("312017" %in% blast_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_complex_summary_curr$arabidopsis <- as.numeric("3702" %in% blast_out_complex_curr$target_taxid) * 0.5
  blast_out_complex_summary_curr$arabidopsis <- blast_out_complex_summary_curr$arabidopsis + as.numeric("3702" %in% blast_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  blast_out_complex_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% blast_out_complex_curr$target_taxid) * 0.5
  blast_out_complex_summary_curr$spinach_chloroplast <- blast_out_complex_summary_curr$spinach_chloroplast + as.numeric("3562" %in% blast_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  
  blast_complex_species <- rbind(blast_complex_species, blast_out_complex_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(blast_complex_species[,-1])

# Binary color
color_mat[color_mat == 0.5] <- 0
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- blast_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$Name))

# Plotting the tile heatmap with ggplot2
pdf("complexI_blast.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# plmsearch analysis
plmsearch_out <- plmsearch_out_BBH
title_string <- "plmsearch BBH similarity score > 0.3 with human homolog"

# Replace swissprot entries with uniprot accessions
swissprot_mapping_subset <- swissprot_mapping[match(plmsearch_out$query_uniprot_entry, swissprot_mapping$uniprot_entry),]
plmsearch_out$query_uniprot_accession <- swissprot_mapping_subset$uniprot_accession
# Convert accessions to uniprot
plmsearch_out$target_uniprot_accession <- orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession[match(plmsearch_out$target, orthogroups_collapse_accessions_seprows_uniprot$accessions)]
# Keep mtDNA accessions (starting with AAK) for Tetrahymena
plmsearch_out$target_uniprot_accession[plmsearch_out$target_taxid == 312017 & grepl("312017_AAK", plmsearch_out$sseqid, fixed=TRUE)] <- plmsearch_out$target_uniprot_entry[plmsearch_out$target_taxid == 312017 & grepl("312017_AAK", plmsearch_out$sseqid, fixed=TRUE)]

plmsearch_out_complex <- plmsearch_out[which(plmsearch_out$query_uniprot_accession %in% complex_table_core$Homo_sapiens),]
plmsearch_out_complex$complex_id <- complex_table_core$Name[match(plmsearch_out_complex$query_uniprot_accession, complex_table_core$Homo_sapiens)]
plmsearch_out_complex <- plmsearch_out_complex %>% rowwise() %>% mutate(true_top_hit = target_uniprot_accession %in% unlist(strsplit(unlist(complex_table_core[which(complex_table_core$Name == complex_id), c("Homo_sapiens", "Tetrahymena_thermophila", "Arabidopsis_thaliana", "Spinach_chloroplast")]), ",")))
plmsearch_out_complex_summary <- plmsearch_out_complex %>% group_by(complex_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))
# plmsearch_out_complex_summary$has_hit <- as.numeric(plmsearch_out_complex_summary$n_hits > 0)

# Mark all Yeast and Giardia hits as true since has been lost
plmsearch_out_complex_summary$any_true_top_hit[plmsearch_out_complex_summary$target_taxid == "559292"] <- TRUE 
plmsearch_out_complex_summary$any_true_top_hit[plmsearch_out_complex_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
plmsearch_complex_species <- c()
for (complex_id_curr in unique(plmsearch_out_complex_summary$complex_id)) {
  plmsearch_out_complex_curr <- plmsearch_out_complex_summary[plmsearch_out_complex_summary$complex_id == complex_id_curr,]
  
  # To use plmsearch true hits
  plmsearch_out_complex_summary_curr_true_hit <- plmsearch_out_complex_curr[plmsearch_out_complex_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  plmsearch_out_complex_summary_curr <- data.frame(complex_id_curr, human = 0, yeast = 0, tetrahymena = 0, arabidopsis = 0, giardia = 0, spinach_chloroplast = 0)
  
  # Binary
  plmsearch_out_complex_summary_curr$human <- as.numeric("9606" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex_summary_curr$yeast <- as.numeric("559292" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex_summary_curr$tetrahymena <- as.numeric("312017" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex_summary_curr$arabidopsis <- as.numeric("3702" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex_summary_curr$giardia <- as.numeric("5741" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid)
  plmsearch_out_complex_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid)
  
  # Tri color
  plmsearch_out_complex_summary_curr$human <- as.numeric("9606" %in% plmsearch_out_complex_curr$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$human <- plmsearch_out_complex_summary_curr$human + as.numeric("9606" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$tetrahymena <- as.numeric("312017" %in% plmsearch_out_complex_curr$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$tetrahymena <- plmsearch_out_complex_summary_curr$tetrahymena + as.numeric("312017" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$arabidopsis <- as.numeric("3702" %in% plmsearch_out_complex_curr$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$arabidopsis <- plmsearch_out_complex_summary_curr$arabidopsis + as.numeric("3702" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% plmsearch_out_complex_curr$target_taxid) * 0.5
  plmsearch_out_complex_summary_curr$spinach_chloroplast <- plmsearch_out_complex_summary_curr$spinach_chloroplast + as.numeric("3562" %in% plmsearch_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  
  plmsearch_complex_species <- rbind(plmsearch_complex_species, plmsearch_out_complex_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(plmsearch_complex_species[,-1])

# Binary color
color_mat[color_mat == 0.5] <- 0
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- plmsearch_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$Name))

# Plotting the tile heatmap with ggplot2
pdf("complexI_plmsearch.plmalign.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()




# foldseek analysis
foldseek_out <- foldseek_out_BBH
title_string <- "foldseek BBH expect<1 with human homolog"

# Convert mtDNA accessions (starting with AAK) for Tetrahymena
tetrahymena_mtdna_mapping <- orthogroups_collapse_accessions_seprows_uniprot[grep("^312017_AAK", orthogroups_collapse_accessions_seprows_uniprot$accessions),]
foldseek_out$target[foldseek_out$target_taxid == 312017 & foldseek_out$target %in% tetrahymena_mtdna_mapping$uniprot_accession] <- gsub("312017_", "", tetrahymena_mtdna_mapping$accessions[match(foldseek_out$target[foldseek_out$target_taxid == 312017 & foldseek_out$target %in% tetrahymena_mtdna_mapping$uniprot_accession], tetrahymena_mtdna_mapping$uniprot_accession)])

foldseek_out_complex <- foldseek_out[which(foldseek_out$query %in% complex_table_core$Homo_sapiens),]
foldseek_out_complex$complex_id <- complex_table_core$Name[match(foldseek_out_complex$query, complex_table_core$Homo_sapiens)]
foldseek_out_complex <- foldseek_out_complex %>% rowwise() %>% mutate(true_top_hit = target %in% c(unlist(strsplit(unlist(complex_table_core[which(complex_table_core$Name == complex_id), c("Homo_sapiens", "Tetrahymena_thermophila", "Arabidopsis_thaliana", "Spinach_chloroplast")]), ",")), ath_idmapping$Entry[which(ath_idmapping$From %in% ath_idmapping$From[match(unlist(strsplit(complex_table_core$Arabidopsis_thaliana[which(complex_table_core$Name == complex_id)], ",")), ath_idmapping$Entry)])]))
foldseek_out_complex_summary <- foldseek_out_complex %>% group_by(complex_id, target_taxid) %>% summarize(n_hits = n(), any_true_top_hit = any(true_top_hit))

# Mark all Yeast and Giardia hits as true since has been lost
foldseek_out_complex_summary$any_true_top_hit[foldseek_out_complex_summary$target_taxid == "559292"] <- TRUE 
foldseek_out_complex_summary$any_true_top_hit[foldseek_out_complex_summary$target_taxid == "5741"] <- TRUE 

# Split into 1 column per species
foldseek_complex_species <- c()
for (complex_id_curr in unique(foldseek_out_complex_summary$complex_id)) {
  foldseek_out_complex_curr <- foldseek_out_complex_summary[foldseek_out_complex_summary$complex_id == complex_id_curr,]
  
  # To use foldseek true hits
  foldseek_out_complex_summary_curr_true_hit <- foldseek_out_complex_curr[foldseek_out_complex_curr$any_true_top_hit,]
  
  # Fractional score for true hit vs any hit
  foldseek_out_complex_summary_curr <- data.frame(complex_id_curr, human = 0, yeast = 0, tetrahymena = 0, arabidopsis = 0, giardia = 0, spinach_chloroplast = 0)
  
  # Binary
  foldseek_out_complex_summary_curr$human <- as.numeric("9606" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid)
  foldseek_out_complex_summary_curr$yeast <- as.numeric("559292" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid)
  foldseek_out_complex_summary_curr$tetrahymena <- as.numeric("312017" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid)
  foldseek_out_complex_summary_curr$arabidopsis <- as.numeric("3702" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid)
  foldseek_out_complex_summary_curr$giardia <- as.numeric("5741" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid)
  foldseek_out_complex_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid)
  
  # Tri color
  foldseek_out_complex_summary_curr$human <- as.numeric("9606" %in% foldseek_out_complex_curr$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$human <- foldseek_out_complex_summary_curr$human + as.numeric("9606" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$tetrahymena <- as.numeric("312017" %in% foldseek_out_complex_curr$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$tetrahymena <- foldseek_out_complex_summary_curr$tetrahymena + as.numeric("312017" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$arabidopsis <- as.numeric("3702" %in% foldseek_out_complex_curr$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$arabidopsis <- foldseek_out_complex_summary_curr$arabidopsis + as.numeric("3702" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$spinach_chloroplast <- as.numeric("3562" %in% foldseek_out_complex_curr$target_taxid) * 0.5
  foldseek_out_complex_summary_curr$spinach_chloroplast <- foldseek_out_complex_summary_curr$spinach_chloroplast + as.numeric("3562" %in% foldseek_out_complex_summary_curr_true_hit$target_taxid) * 0.5
  
  foldseek_complex_species <- rbind(foldseek_complex_species, foldseek_out_complex_summary_curr)
}

# Reshape the color data into a dataframe for plotting
color_mat <- as.matrix(foldseek_complex_species[,-1])

# Binary color
color_mat[color_mat == 0.5] <- 0
color_mat <- ifelse(color_mat, "#0047AB", "white")

color_df <- as.data.frame(color_mat)

color_df$row <- foldseek_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$Name))

# Plotting the tile heatmap with ggplot2
pdf("complexI_foldseek.BBH_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = title_string) +
  theme_minimal()
dev.off()



## EggNOG complex I
# Idmapping for human and yeast complex
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "complexI_human_uniprot_to_ensemblprotein_idmapping_2025_07_10.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

eggnog_species_ids <- c("9606", "312017", "3702")

OG_complex <- c()
for (i in 1:nrow(complex_table_core)) {
  complex_curr_row <- complex_table_core[i,]
  
  # Get all Arabidopsis paralogs
  ath_ids <- unlist(strsplit(complex_curr_row$Arabidopsis_thaliana, ","))
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid_curr <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- paste0("3702.", ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid_curr)])
  
  # Map ids to eggnog
  ath_geneid <- paste0("3702.", ath_ids)
  tth_geneid <- paste0("312017.", complex_curr_row$Tetrahymena_thermophila)
  spinacea_geneid <- paste0("3562.", complex_curr_row$Spinach_chloroplast)
  hsa_geneid <- paste0("9606.", idmapping_human$To[idmapping_human$From %in% complex_curr_row$Homo_sapiens])
  
  complex_curr_accessions <- c(hsa_geneid, sce_geneid, ath_geneid, tth_geneid, spinacea_geneid, ath_uniprot_ids)
  
  OG_curr_row <- root_OG_long %>% filter(protein_ids %in% complex_curr_accessions)
  OG_curr_row$complex_id <- complex_curr_row$Name
  
  # Include homologs from other species
  OG_curr_all_members <- root_OG_long[root_OG_long$OG_id %in% OG_curr_row$OG_id,]
  yeast_giardia_index <- which(gsub("\\..*", "", OG_curr_all_members$protein_ids) %in% c("4932", "184922"))
  if (length(yeast_giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[yeast_giardia_index,]
    OG_curr_all_members$complex_id <- complex_curr_row$Name
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "protein_ids", "complex_id")]
  
  # Find missing entries
  missing_species <- eggnog_species_ids[!eggnog_species_ids %in% gsub("\\..*", "", OG_curr_row$protein_ids)]
  if (length(missing_species) > 0) {
    print(complex_curr_row$Name)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", protein_ids = missing_species, complex_id = complex_curr_row$Name)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_complex <- rbind(OG_complex, OG_curr_row)
}

# Get taxids
OG_complex$taxid <- gsub("\\..*", "", OG_complex$protein_ids)

OG_complex_summary <- OG_complex %>% group_by(complex_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_complex_species <- c()
for (complex_id_curr in unique(OG_complex_summary$complex_id)) {
  OG_complex_summary_curr <- OG_complex_summary[OG_complex_summary$complex_id == complex_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 4932]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  tetrahymena_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 312017]
  if (length(tetrahymena_OG_id) == 0) {
    tetrahymena_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  giardia_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 184922]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  spinach_chloroplast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_OG_id) == 0) {
    spinach_chloroplast_OG_id <- NA
  }
  
  OG_complex_species_curr <- data.frame(complex_id_curr, `9606` = human_OG_id, `559292` = yeast_OG_id, `185431` = tetrahymena_OG_id, `3702` = arabidopsis_OG_id, `5741` = giardia_OG_id, `3562` = spinach_chloroplast_OG_id, check.names = FALSE)
  OG_complex_species <- rbind(OG_complex_species, OG_complex_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_complex_species_bool <- t(apply(OG_complex_species[,-1], 1, function(row) row == row[1]))
OG_complex_species_bool[is.na(OG_complex_species_bool)] <- FALSE
OG_complex_species_colors <- ifelse(OG_complex_species_bool, "#0047AB", "white")
OG_complex_species_colors <- as.data.frame(OG_complex_species_colors)

# Mark proteins missing from dataset as gray
OG_complex_species_colors[which(OG_complex_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_complex_species_colors

color_df$row <- OG_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$Name))


# Plotting the tile heatmap with ggplot2
pdf("complexI_EggNOG_Eukaryota_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "EggNOG") +
  theme_minimal()
dev.off()



## OrthoDB complex I

# Idmapping for human and yeast complex
idmapping_human <- read.table(here("data", "annotation", "pathways", "idmapping", "complexI_human_uniprot_to_ensembl_idmapping_2025_07_10.tsv"), sep="\t", header=TRUE)
idmapping_human$To <- gsub("\\..*", "", idmapping_human$To)

orthodb_species_ids <- c("312017_0", "3702_0", "9606_0")

OG_complex <- c()
for (i in 1:nrow(complex_table_core)) {
  complex_curr_row <- complex_table_core[i,]
  
  hsa_ensg_id <- idmapping_human$To[idmapping_human$From %in% complex_curr_row$Homo_sapiens]
  
  # Get all Arabidopsis paralogs
  ath_ids <- unlist(strsplit(complex_curr_row$Arabidopsis_thaliana, ","))
  
  # Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
  ath_geneid <- ath_idmapping$From[match(ath_ids, ath_idmapping$Entry)]
  ath_uniprot_ids <- ath_idmapping$Entry[which(ath_idmapping$From == ath_geneid)]
  idmapping_ath_curr <- idmapping_ath %>% filter(KEGG %in% ath_geneid | TAIR %in% ath_geneid | Araport %in% ath_geneid)
  ath_uniprot_ids <- c(ath_uniprot_ids, idmapping_ath_curr$From)
  
  complex_curr_accessions <- c(complex_curr_row$Tetrahymena_thermophila, ath_ids, complex_curr_row$Spinach_chloroplast, ath_uniprot_ids)
  complex_curr_accessions <- complex_curr_accessions[complex_curr_accessions != ""]
  
  OG_curr_row <- odb_eukaryota_ogs_genes %>% filter(uniprot_id %in% complex_curr_accessions | gene_id %in% hsa_ensg_id | gene_id %in% ath_geneid)
  OG_curr_row$complex_id <- complex_curr_row$Name
  
  # Include homologs from other species
  OG_curr_all_members <- odb_eukaryota_ogs_genes[odb_eukaryota_ogs_genes$OG_id %in% OG_curr_row$OG_id,]
  yeast_giardia_index <- which(OG_curr_all_members$species_id %in% c("559292_0", "5741_0"))
  if (length(yeast_giardia_index) > 0) {
    OG_curr_all_members <- OG_curr_all_members[yeast_giardia_index,]
    OG_curr_all_members$complex_id <- complex_curr_row$Name
    OG_curr_row <- rbind(OG_curr_row, OG_curr_all_members)
  }
  
  OG_curr_row <- OG_curr_row[,c("OG_id", "species_id", "complex_id")]
  
  # Find missing entries
  missing_species <- orthodb_species_ids[!orthodb_species_ids %in% OG_curr_row$species_id]
  if (length(missing_species) > 0) {
    print(complex_curr_row$Name)
    print(missing_species)
    missing_rows <- data.frame(OG_id = "protein_missing_from_dataset", species_id = missing_species, complex_id = complex_curr_row$Name)
    OG_curr_row <- rbind(OG_curr_row, missing_rows)
  }
  
  OG_complex <- rbind(OG_complex, OG_curr_row)
}

# Get taxids
OG_complex$taxid <- gsub("_.*", "", OG_complex$species_id)

OG_complex_summary <- OG_complex %>% group_by(complex_id, taxid) %>% summarize(OG_id = paste0(unique(OG_id), collapse=","))

# Split into 1 column per species
OG_complex_species <- c()
for (complex_id_curr in unique(OG_complex_summary$complex_id)) {
  OG_complex_summary_curr <- OG_complex_summary[OG_complex_summary$complex_id == complex_id_curr,]
  
  # Get OG ids
  human_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 9606]
  if (length(human_OG_id) == 0) {
    human_OG_id <- NA
  }
  yeast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 559292]
  if (length(yeast_OG_id) == 0) {
    yeast_OG_id <- NA
  }
  tetrahymena_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 312017]
  if (length(tetrahymena_OG_id) == 0) {
    tetrahymena_OG_id <- NA
  }
  arabidopsis_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3702]
  if (length(arabidopsis_OG_id) == 0) {
    arabidopsis_OG_id <- NA
  }
  giardia_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 5741]
  if (length(giardia_OG_id) == 0) {
    giardia_OG_id <- NA
  }
  spinach_chloroplast_OG_id <- OG_complex_summary_curr$OG_id[OG_complex_summary_curr$taxid == 3562]
  if (length(spinach_chloroplast_OG_id) == 0) {
    spinach_chloroplast_OG_id <- NA
  }
  
  OG_complex_species_curr <- data.frame(complex_id_curr, `9606` = human_OG_id, `559292` = yeast_OG_id, `312017` = tetrahymena_OG_id, `3702` = arabidopsis_OG_id, `5741` = giardia_OG_id, `3562` = spinach_chloroplast_OG_id, check.names = FALSE)
  OG_complex_species <- rbind(OG_complex_species, OG_complex_species_curr)
}

# Plot heatmap
# Function to assign colors to each row
assign_colors <- function(row) {
  # Initialize the color vector with the first color as blue
  colors <- c("blue")
  
  # Map each unique value in the row to a color
  color_map <- setNames("blue", row[1]) # Start with the first value being blue
  color_list <- c("red", "purple", "darkgreen") # other colors to assign
  color_index <- 1
  
  for (i in 2:length(row)) {
    value <- row[i]
    # If the value is already mapped, use the same color
    if (value %in% names(color_map)) {
      colors <- c(colors, color_map[value])
    } else {
      # Assign a new color to this unique value
      color_map[value] <- color_list[color_index]
      colors <- c(colors, color_list[color_index])
      color_index <- color_index + 1
    }
  }
  
  return(colors)
}

# Only blue/white
OG_complex_species_bool <- t(apply(OG_complex_species[,-1], 1, function(row) row == row[1]))
OG_complex_species_bool[is.na(OG_complex_species_bool)] <- FALSE
OG_complex_species_colors <- ifelse(OG_complex_species_bool, "#0047AB", "white")
OG_complex_species_colors <- as.data.frame(OG_complex_species_colors)

# Mark proteins missing from dataset as gray
OG_complex_species_colors[which(OG_complex_species[,-1] == "protein_missing_from_dataset", arr.ind = TRUE)] <- "gray"

# Reshape the color data into a dataframe for plotting
color_df <- OG_complex_species_colors

color_df$row <- OG_complex_species$complex_id_curr  # Add row identifier

# Melt data for ggplot2
color_long <- melt(color_df, id.vars = "row", variable.name = "column", value.name = "color")
color_long$row <- factor(color_long$row, levels = rev(complex_table_core$Name))


# Plotting the tile heatmap with ggplot2
pdf("complexI_OrthoDB_heatmap.pdf", width = 6, height = 8)
ggplot(color_long, aes(x = column, y = row)) +
  geom_tile(aes(fill = color), color = "white") +
  scale_fill_identity() +  # Use colors directly from the "color" column
  labs(x = "", y = "", title = "OrthoDB") +
  theme_minimal()
dev.off()
