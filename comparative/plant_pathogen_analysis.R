
library(tidyverse)
library(ape)
library(phytools)
library(TreeTools)

library(ggplot2)
library(ggrepel)
library(gplots)
library(Matrix)

# Set ggplot theme
theme_set(theme_classic())

# Get absense predictions for Arabidopsis
homology_power_per_species <- read.table(here("data", "abSENSE_HMM", "species_tree_1", "absense_results_per_species.tsv"), sep="\t", header=FALSE)
colnames(homology_power_per_species) <- c("PhROG_id", "taxid", "probability_of_detection")
homology_power_per_species <- homology_power_per_species %>% rowwise() %>% mutate(OG_id = gsub("_Node.*", "", PhROG_id))
homology_power_ath <- homology_power_per_species %>% filter(taxid == "3702")

gold_gene_accession_OG_id_df_primary <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2026.04.05.tsv"), sep="\t", header=TRUE)
primary_mito_OGs_mitocarta_or_mtdna <- unique(unlist(strsplit(gold_gene_accession_OG_id_df_primary$OG_id, split=",")))

puccinia_taxid <- "56615"
fusarium_taxid <- "660029"
colletotrichum_taxid <- "1209926"
fungal_pathogens <- c(puccinia_taxid, fusarium_taxid, colletotrichum_taxid)

ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk673spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")
ogs_long_primary <- ogs_long %>% filter(BOOL_PRIMARY_OG)

uniprot_proteomes_203euks_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

# Read in DeepLoc predictions
deeploc_results <- read.table(here("data/deeploc/predictions", "DeepLoc2.0-mito_predictions.tsv"), header=TRUE)
colnames(deeploc_results) <- c("Protein_ID", "Mitochondrion")
deeploc_results$taxid <- gsub("_.*", "", deeploc_results$Protein_ID)

n_mito_species_threshold <- 5

deeploc_thresholds <- read.csv(here("data/deeploc", "deeploc_thresholds.csv"))
deeploc_mito_threshold <- deeploc_thresholds$threshold[deeploc_thresholds$label == "Mitochondrion"]
mito_threshold <- deeploc_mito_threshold

deeploc_results$label <- "Nonmito"
deeploc_results$label[deeploc_results$Mitochondrion >= mito_threshold] <- "Mito"

deeploc_results <- merge(deeploc_results, ogs_long, by.x=c("Protein_ID", "taxid"), by.y=c("accession", "taxid"), all.x=TRUE)

# Mark organelle-encoded proteins, if available
all_nonmito_accessions <- read.table(here("data/deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
all_mito_accessions <- read.table(here("data/deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1
deeploc_results$label[deeploc_results$Protein_ID %in% all_nonmito_accessions] <- "Nonmito"
deeploc_results$label[deeploc_results$Protein_ID %in% all_mito_accessions] <- "Mito"

# Ignore species missing organelle genomes in OGs that have organelle encoded proteins
missing_nonmito_organelle_taxids <- read.table(here("data/deeploc", "missing_nonmito_organelle_taxids.txt"))$V1
missing_mtdna_taxids <- read.table(here("data/deeploc", "missing_mtdna_taxids.txt"))$V1
missing_mtdna_taxids <- c(missing_mtdna_taxids, puccinia_taxid, colletotrichum_taxid)
nonmito_organelle_OG_ids <- unique(deeploc_results$Orthogroup[deeploc_results$Protein_ID %in% all_nonmito_accessions])
deeploc_results$label[deeploc_results$Orthogroup %in% nonmito_organelle_OG_ids & deeploc_results$taxid %in% missing_nonmito_organelle_taxids] <- NA
mito_organelle_OG_ids <- unique(deeploc_results$Orthogroup[deeploc_results$Protein_ID %in% all_mito_accessions])
deeploc_results$label[deeploc_results$Orthogroup %in% mito_organelle_OG_ids & deeploc_results$taxid %in% missing_mtdna_taxids] <- NA

gold_gene_accession_OG_id_df <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_2026.04.05.tsv"), sep="\t", header=TRUE)
gold_gene_accession_OG_id_df$taxid <- gsub("_.*", "", gold_gene_accession_OG_id_df$gene_accession)

deeploc_results$type <- "Predicted"
complete_mitoproteome_species <- c("9606", "559292", "3702", "1257118", "5741", "5689", "185431", "32595")
deeploc_results$type[deeploc_results$taxid %in% complete_mitoproteome_species] <- "Experimental"
deeploc_results$label[deeploc_results$taxid %in% complete_mitoproteome_species] <- "Nonmito"
deeploc_results$label[deeploc_results$Protein_ID %in% gold_gene_accession_OG_id_df$gene_accession[gold_gene_accession_OG_id_df$taxid %in% complete_mitoproteome_species]] <- "Mito"

## Find OGs that have at least 5 species with predicted mito protein. Using primary OGs
deeploc_results_gte5spp <- deeploc_results %>% filter(taxid %in% uniprot_proteomes_203euks_tax$TaxId) %>% group_by(Orthogroup) %>% filter(label == "Mito") %>% filter(length(unique(taxid)) >= n_mito_species_threshold)
deeploc_results_gte5spp_OG_ids <- unique(c(deeploc_results_gte5spp$Orthogroup, primary_mito_OGs_mitocarta_or_mtdna))

# Mark low confidence predictions as nonmito
deeploc_results$label[!deeploc_results$Orthogroup %in% deeploc_results_gte5spp_OG_ids & deeploc_results$label == "Mito"] <- "Nonmito"

deeploc_results_primary <- deeploc_results %>% filter(BOOL_PRIMARY_OG)

## Plant parasite targets via OGs
puccinia_taxid <- "56615"
fusarium_taxid <- "660029"
colletotrichum_taxid <- "1209926"
ustilago_taxid <- "237631"
fungal_pathogen_taxids <- c(puccinia_taxid, fusarium_taxid, colletotrichum_taxid, ustilago_taxid)
phytophthora_infestans_taxid <- "4787"
plasmodiophora_taxid <- "37360"
plant_pathogen_taxids <- c(fungal_pathogen_taxids, phytophthora_infestans_taxid, plasmodiophora_taxid)
selected_taxids <- c("3702", phytophthora_infestans_taxid, plasmodiophora_taxid, puccinia_taxid, ustilago_taxid, colletotrichum_taxid, fusarium_taxid, "559292", "9606")

ogs_long_plant_pathogen <- ogs_long %>% filter(taxid %in% selected_taxids) %>% group_by(Orthogroup) %>% mutate(fungi = as.numeric(any(fungal_pathogen_taxids %in% taxid)), phytophthora = as.numeric(phytophthora_infestans_taxid %in% taxid), plasmodiophora = as.numeric(plasmodiophora_taxid %in% taxid)) %>% mutate(present_n_parasite_clades = fungi + phytophthora + plasmodiophora) %>% filter(present_n_parasite_clades >= 2)

# Assign mito labels
ogs_long_plant_pathogen$Mitochondrion <- deeploc_results_primary$Mitochondrion[match(ogs_long_plant_pathogen$accession, deeploc_results_primary$Protein_ID)]
ogs_long_plant_pathogen$label <- deeploc_results_primary$label[match(ogs_long_plant_pathogen$accession, deeploc_results_primary$Protein_ID)]

ogs_long_plant_pathogen <- ogs_long_plant_pathogen %>% filter(label == "Mito") %>% group_by(Orthogroup) %>% mutate(fungi = as.numeric(any(fungal_pathogen_taxids %in% taxid)), phytophthora = as.numeric(phytophthora_infestans_taxid %in% taxid), plasmodiophora = as.numeric(plasmodiophora_taxid %in% taxid)) %>% mutate(present_n_parasite_clades = fungi + phytophthora + plasmodiophora) %>% filter(present_n_parasite_clades >= 1)

ogs_long_ath <- ogs_long %>% filter(taxid == "3702")
ath_OG_ids <- unique(ogs_long_ath$Orthogroup)
ogs_long_plant_pathogen$present_ath_nonmito <- as.numeric(ogs_long_plant_pathogen$Orthogroup %in% ath_OG_ids)

plant_pathogen_OG_ids <- unique(ogs_long_plant_pathogen$Orthogroup)


# Mark homologs as probability of detection = 1
ogs_long_select <- ogs_long %>% filter(taxid %in% "3702") %>% filter(Orthogroup %in% plant_pathogen_OG_ids)
homology_power_ath <- homology_power_ath %>% mutate(bool_has_homolog = taxid %in% ogs_long_select$taxid[ogs_long_select$Orthogroup == OG_id])
homology_power_ath_detectable <- homology_power_ath %>% filter(probability_of_detection > 0.95 | bool_has_homolog)
plant_pathogen_OG_ids_detectable_ath <- unique(homology_power_ath_detectable$OG_id)
ogs_long_plant_pathogen$absense <- "Not powered"
ogs_long_plant_pathogen$absense[ogs_long_plant_pathogen$Orthogroup %in% plant_pathogen_OG_ids_detectable_ath] <- "Powered"


# Get PFAM annotations
pfam_seprows <- read.table(here("data", "annotation", "pfam", "euk673spp_prokgroups_pfam_list.filtered_seprows.tsv"), sep="\t", header=TRUE)
ogs_long_plant_pathogen_pfam <- merge(ogs_long_plant_pathogen, pfam_seprows, by.x = "accession", by.y = "Accession", all.x = TRUE)
ogs_long_plant_pathogen_pfam <- ogs_long_plant_pathogen_pfam %>% group_by(accession, Orthogroup, taxid, BOOL_PRIMARY_OG, fungi, phytophthora, plasmodiophora, present_n_parasite_clades, Mitochondrion, label, present_ath_nonmito, absense) %>% summarize(PfamList = paste0(PfamList, collapse=","), .groups = "drop")
ogs_long_plant_pathogen_pfam$PfamList[is.na(ogs_long_plant_pathogen_pfam$PfamList)] <- ""
ogs_long_plant_pathogen_pfam <- ogs_long_plant_pathogen_pfam[order(ogs_long_plant_pathogen_pfam$Orthogroup, ogs_long_plant_pathogen_pfam$taxid),]

ogs_long_plant_pathogen_pfam_no_ath <- ogs_long_plant_pathogen_pfam %>% filter(present_ath_nonmito == 0)
length(unique(ogs_long_plant_pathogen_pfam_no_ath$Orthogroup))

ogs_long_plant_pathogen_pfam_no_ath_and_detectable <- ogs_long_plant_pathogen_pfam_no_ath %>% filter(absense == "Powered")
length(unique(ogs_long_plant_pathogen_pfam_no_ath_and_detectable$Orthogroup))

ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared <- ogs_long_plant_pathogen_pfam_no_ath_and_detectable %>% filter(present_n_parasite_clades >= 2)
length(unique(ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared$Orthogroup))

ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared <- ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared %>% filter(taxid %in% selected_taxids)


# Venn diagram on a per PhROG basis
parasite_mito_set <- list('Fungi' = unique(ogs_long_plant_pathogen$Orthogroup[ogs_long_plant_pathogen$fungi == 1]), 'Phytophthora' = unique(ogs_long_plant_pathogen$Orthogroup[ogs_long_plant_pathogen$phytophthora == 1]), 'Plasmodiophora' = unique(ogs_long_plant_pathogen$Orthogroup[ogs_long_plant_pathogen$plasmodiophora == 1]))
library(ggvenn)
pdf('phytopathogen_PhROGs_venn.pdf', width=6, height=4)
ggvenn(parasite_mito_set, c("Fungi", "Phytophthora", "Plasmodiophora"), show_percentage = FALSE, set_name_size=4, text_size=6, fill_color = c("gray", "#00A36C", "purple"))
dev.off()

parasite_noath_powered_mito_set <- list('Fungi' = unique(ogs_long_plant_pathogen$Orthogroup[ogs_long_plant_pathogen$fungi == 1 & ogs_long_plant_pathogen$present_ath_nonmito == 0 & ogs_long_plant_pathogen$absense == "Powered"]), 'Phytophthora' = unique(ogs_long_plant_pathogen$Orthogroup[ogs_long_plant_pathogen$phytophthora == 1 & ogs_long_plant_pathogen$present_ath_nonmito == 0 & ogs_long_plant_pathogen$absense == "Powered"]), 'Plasmodiophora' = unique(ogs_long_plant_pathogen$Orthogroup[ogs_long_plant_pathogen$plasmodiophora == 1 & ogs_long_plant_pathogen$present_ath_nonmito == 0 & ogs_long_plant_pathogen$absense == "Powered"]))
pdf('phytopathogen_PhROGs_no.Arabidopsis.powered_venn.pdf', width=6, height=4)
ggvenn(parasite_noath_powered_mito_set, c("Fungi", "Phytophthora", "Plasmodiophora"), show_percentage = FALSE, set_name_size=4, text_size=6, fill_color = c("gray", "#00A36C", "purple"))
dev.off()


## Plot heatmap
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

# Build OGs matrix
OG_ids <- unique(ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared$Orthogroup)
ogs_long_subset <- ogs_long %>% filter(taxid %in% selected_taxids)
og_presence_mat <- matrix(data = 0, nrow = length(OG_ids), ncol = length(selected_taxids))
colnames(og_presence_mat) <- selected_taxids
symbols <- c()
for (i in 1:length(OG_ids)) {
  curr_OG_id <- OG_ids[i]
  ogs_subset_long_curr <- ogs_long_subset %>% filter(Orthogroup == curr_OG_id)
  ogs_subset_long_curr_nonmito <- ogs_subset_long_curr %>% filter(!accession %in% ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared$accession)
  ogs_subset_long_curr_mito <- ogs_subset_long_curr %>% filter(accession %in% ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared$accession)
  og_presence_mat[i, unique(ogs_subset_long_curr_nonmito$taxid)] <- 0.5
  og_presence_mat[i, unique(ogs_subset_long_curr_mito$taxid)] <- 1
  
  unique_pfams <- unique(unlist(strsplit(ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared$PfamList[ogs_long_plant_pathogen_pfam_no_ath_and_detectable_shared$Orthogroup %in% curr_OG_id], ",")))
  unique_pfams <- unique_pfams[unique_pfams != "NA"]
  curr_symbol <- paste0(gsub("_.*", "", curr_OG_id), ":", paste0(sort(unique_pfams), collapse=","))
  if (curr_symbol != "") {
    symbols <- c(symbols, curr_symbol)
  } else {
    symbols <- c(symbols, curr_OG_id)
  }
}
rownames(og_presence_mat) <- symbols

og_presence_mat <- og_presence_mat[,selected_taxids]

# Sort matrix for plotting
og_presence_mat_plot <- og_presence_mat
rownames(og_presence_mat_plot) <- str_trunc(rownames(og_presence_mat), 40)

# Order by values
n_parasite_kingdoms <- as.numeric(rowSums(og_presence_mat_plot[,fungal_pathogen_taxids]) == 4) + as.numeric(og_presence_mat_plot[,phytophthora_infestans_taxid] == 1) + as.numeric(og_presence_mat_plot[,plasmodiophora_taxid] == 1)
og_presence_mat_reorder <- og_presence_mat_plot[order(og_presence_mat_plot[,"3702"], -n_parasite_kingdoms, -rowSums(og_presence_mat_plot)),]

colnames(og_presence_mat_reorder) <- uniprot_proteomes_tax$ScientificName[match(colnames(og_presence_mat_reorder), uniprot_proteomes_tax$TaxId)]

my_palette <- colorRampPalette(c("white", "#999999", "#eb35ff"))(n = 100)
pdf(paste0("phytopathogen_targets_heatmap.pdf"), height=10, width=8)
heatmap.2(og_presence_mat_reorder, Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none", margins=c(16,16), na.rm=TRUE, col=my_palette, trace="none", key.xlab="Absence/Presence", cexRow = 1, cexCol = 1, key=FALSE, offsetRow=0)
dev.off()

