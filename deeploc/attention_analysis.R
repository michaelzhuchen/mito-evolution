
library(tidyverse)
library(ggplot2)
library(reshape2)
theme_set(theme_classic())


# Read in datasets
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")
ogs_long_primary <- ogs_long %>% filter(BOOL_PRIMARY_OG)
ogs_long_primary_bacteria <- ogs_long_primary %>% filter(taxid %in% uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Bacteria"])
ogs_long_primary_archaea <- ogs_long_primary %>% filter(taxid %in% uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Archaea"])

origin_table <- read.table(here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t")
colnames(origin_table) <- c("OG_id", "n_euk_species_largest_euk_clade", "n_prok_species_largest_prok_clade", "LCA_node", "LCA_node_euks", "origin_domain", "dropped_tips")
origin_table$OG_id <- gsub(".faa_clipkit.gappy.msa", "", origin_table$OG_id, fixed=TRUE)
prok_origin_OG_ids <- origin_table$OG_id[origin_table$origin_domain == "Prokaryote"]
euk_origin_OG_ids <- origin_table$OG_id[origin_table$origin_domain == "Eukaryote"]

gold_gene_accession_OG_id_df <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2026.04.05.tsv"), sep="\t", header=TRUE)

all_nonmito_accessions <- read.table(here("data/deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
all_mtdna_accessions <- read.table(here("data/deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1

# Get TargetP results
targetp_results <- read.table(here("data", "annotation", "targetp", "euk203spp.mtDNA_prokgroups_underscore_targetp2.0.tsv"), sep="\t", header=TRUE)
targetp_results <- targetp_results %>% mutate(taxid = gsub("_.*", "", ID))
targetp_results$Orthogroup <- ogs_long_primary$Orthogroup[match(targetp_results$ID, ogs_long_primary$accession)]
# Remove orthogroups with any mtDNA encoded proteins
targetp_results <- targetp_results %>% filter(!ID %in% c(all_mito_accessions, all_nonmito_accessions))
targetp_results <- targetp_results %>% mutate(has_MTS = (Prediction == "mTP"))

# Get human MitoCarta targetp results
human_mitocarta <- read.csv(here("data", "annotation", "human_mitocarta", "Human.MitoCarta3.0.sheetA.csv"))
human_id2training <- read.table(here("data", "orthogroups", "idmapping", "human.id2training_uniprot_symbols.txt"), sep="\t", header=TRUE)
human_mitocarta <- merge(human_id2training, human_mitocarta, by=c("HumanGeneID"))
human_mitocarta$gene_accession <- paste0("9606_", human_mitocarta$Entry)
human_homologs <- read.delim(here("data", "annotation", "human_mitocarta", "human.homologs.txt"))
human_homologs <- merge(human_id2training, human_homologs, by=c("HumanGeneID"))
human_homologs$gene_accession <- paste0("9606_", human_homologs$Entry)
targetp_results_human_mito <- targetp_results %>% mutate(taxid = gsub("_.*", "", ID)) %>% filter(taxid %in% "9606") %>% filter(ID %in% gold_gene_accession_OG_id_df$gene_accession)

mitofates_threshold <- 0.385
targetp_results_human_mito$MitoFates <- human_homologs$MitoFates1.2_MTS_probability[match(targetp_results_human_mito$ID, human_homologs$gene_accession)]
targetp_results_human_mito$MitoCarta3.0_SubMitoLocalization <- human_mitocarta$MitoCarta3.0_SubMitoLocalization[match(targetp_results_human_mito$ID, human_mitocarta$gene_accession)]

# TargetP only
targetp_results_human_mito <- targetp_results_human_mito %>% mutate(has_MTS = Prediction == "mTP", prokaryote_origin = Orthogroup %in% prok_origin_OG_ids, eukaryote_origin = Orthogroup %in% euk_origin_OG_ids)

selected_taxids <- c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595")
targetp_results_selected <- targetp_results %>% mutate(taxid = gsub("_.*", "", ID)) %>% filter(taxid %in% selected_taxids) %>% filter(ID %in% gold_gene_accession_OG_id_df$gene_accession)


# Compare MTS fraction of extant species
possible_mts_threshold <- 0.2

selected_taxids <- rev(c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595"))
targetp_results_selected <- targetp_results %>% filter(taxid %in% selected_taxids) %>% filter(ID %in% gold_gene_accession_OG_id_df$gene_accession)

# Assign species names
targetp_results_selected$species <- factor(uniprot_proteomes_tax$ScientificName[match(targetp_results_selected$taxid, uniprot_proteomes_tax$TaxId)], levels = uniprot_proteomes_tax$ScientificName[match(rev(selected_taxids), uniprot_proteomes_tax$TaxId)])

selected_columns <- c("ID", "Orthogroup", "taxid", "species", "mTP", "Prediction")
targetp_results_reduce <- rbind(targetp_results_selected[,selected_columns])

# Plot number of MTS/no MTS proteins
targetp_results_reduce$label <- "No MTS"
targetp_results_reduce$label[targetp_results_reduce$mTP > possible_mts_threshold] <- "Possible MTS"
targetp_results_reduce$label[targetp_results_reduce$Prediction == "mTP"] <- "MTS"
targetp_results_reduce$label <- factor(targetp_results_reduce$label, levels = c("No MTS", "Possible MTS", "MTS"))

pdf('barplot_number_MTS_proteins.pdf', height = 6, width = 6)
ggplot(data=targetp_results_reduce, aes(x = species, fill = label)) + geom_bar() + ylab("Mitochondrial proteins") + scale_fill_manual(values = c("MTS" = "black", "Possible MTS" = "#999999", "No MTS" = "#00A36C")) + ylim(0,2000) + coord_flip()
dev.off()

# Plot fraction of MTS/no MTS proteins.
targetp_results_reduce_fraction_MTS <- targetp_results_reduce %>% group_by(species) %>% summarize(fraction = sum(Prediction == "mTP") / n()) %>% mutate(label = "MTS")
targetp_results_reduce_fraction_possibleMTS <- targetp_results_reduce %>% group_by(species) %>% summarize(fraction = sum(Prediction != "mTP" & mTP > possible_mts_threshold) / n()) %>% mutate(label = "Possible MTS")
targetp_results_reduce_fraction_noMTS <- targetp_results_reduce %>% group_by(species) %>% summarize(fraction = sum(!(Prediction == "mTP" | mTP > possible_mts_threshold)) / n()) %>% mutate(label = "No MTS")
targetp_results_reduce_fraction_combined <- rbind(targetp_results_reduce_fraction_MTS, targetp_results_reduce_fraction_possibleMTS, targetp_results_reduce_fraction_noMTS)
targetp_results_reduce_fraction_combined$label <- factor(targetp_results_reduce_fraction_combined$label, levels = c("No MTS", "Possible MTS", "MTS"))
pdf('barplot_fraction_MTS_proteins.pdf', height = 6, width = 6)
ggplot(data=targetp_results_reduce_fraction_combined, aes(x = species, y = fraction, fill = label)) + geom_col() + ylab("Fraction of mitochondrial proteins") + scale_fill_manual(values = c("MTS" = "black", "Possible MTS" = "#999999", "No MTS" = "#00A36C")) + coord_flip()
dev.off()

# Compare MTS fraction of human MitoCarta by subcompartment
selected_compartment <- c("Matrix", "MIM", "IMS", "MOM", "unknown")
targetp_results_human_mito_filter <- targetp_results_human_mito %>% filter(MitoCarta3.0_SubMitoLocalization %in% selected_compartment)
targetp_results_human_mito_filter$MitoCarta3.0_SubMitoLocalization <- factor(targetp_results_human_mito_filter$MitoCarta3.0_SubMitoLocalization, levels = rev(selected_compartment))

targetp_results_human_mito_filter$label <- "No MTS"
targetp_results_human_mito_filter$label[targetp_results_human_mito_filter$mTP > possible_mts_threshold] <- "Possible MTS"
targetp_results_human_mito_filter$label[targetp_results_human_mito_filter$Prediction == "mTP"] <- "MTS"
targetp_results_human_mito_filter$label <- factor(targetp_results_human_mito_filter$label, levels = c("No MTS", "Possible MTS", "MTS"))

# Plot number of MTS/no MTS proteins
pdf('barplot_number_MTS_proteins_human.mitocarta.subcompartments.pdf', height = 6, width = 6)
ggplot(data=targetp_results_human_mito_filter, aes(x = MitoCarta3.0_SubMitoLocalization, fill = label)) + geom_bar() + ylab("Mito proteins") + scale_fill_manual(values = c("MTS" = "black", "Possible MTS" = "#999999", "No MTS" = "#00A36C")) + ylim(0,600) + coord_flip()
dev.off()

# Plot fraction of MTS/no MTS proteins.
targetp_results_human_mito_filter_fraction_MTS <- targetp_results_human_mito_filter %>% group_by(MitoCarta3.0_SubMitoLocalization) %>% summarize(fraction = sum(Prediction == "mTP") / n()) %>% mutate(label = "MTS")
targetp_results_human_mito_filter_fraction_possibleMTS <- targetp_results_human_mito_filter %>% group_by(MitoCarta3.0_SubMitoLocalization) %>% summarize(fraction = sum(Prediction != "mTP" & mTP > possible_mts_threshold) / n()) %>% mutate(label = "Possible MTS")
targetp_results_human_mito_filter_fraction_noMTS <- targetp_results_human_mito_filter %>% group_by(MitoCarta3.0_SubMitoLocalization) %>% summarize(fraction = sum(!(Prediction == "mTP" | mTP > possible_mts_threshold)) / n()) %>% mutate(label = "No MTS")
targetp_results_human_mito_filter_fraction_combined <- rbind(targetp_results_human_mito_filter_fraction_MTS, targetp_results_human_mito_filter_fraction_possibleMTS, targetp_results_human_mito_filter_fraction_noMTS)
targetp_results_human_mito_filter_fraction_combined$label <- factor(targetp_results_human_mito_filter_fraction_combined$label, levels = c("No MTS", "Possible MTS", "MTS"))
pdf('barplot_fraction_MTS_proteins_human.mitocarta.subcompartments.pdf', height = 6, width = 6)
ggplot(data=targetp_results_human_mito_filter_fraction_combined, aes(x = MitoCarta3.0_SubMitoLocalization, y = fraction, fill = label)) + geom_col() + ylab("Fraction of mito. proteome") + scale_fill_manual(values = c("MTS" = "black", "Possible MTS" = "#999999", "No MTS" = "#00A36C")) + coord_flip()
dev.off()



### DeepLoc attention plots
bin_width <- 0.01
selected_taxids <- c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595")
selected_submito_compartments <- c("Matrix", "MIM", "IMS", "MOM")

attn_human <- read.csv(here("data", "deeploc", "attention", "9606.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_yeast <- read.csv(here("data", "deeploc", "attention", "559292.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_aca <- read.csv(here("data", "deeploc", "attention", "1257118.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_gla <- read.csv(here("data", "deeploc", "attention", "5741.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_ath <- read.csv(here("data", "deeploc", "attention", "3702.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_bdi <- read.csv(here("data", "deeploc", "attention", "32595.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_tbr <- read.csv(here("data", "deeploc", "attention", "185431.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)
attn_lta <- read.csv(here("data", "deeploc", "attention", "5689.fasta.split_retrained_deeploc_All_Accurate_attention_combined.csv"), header=FALSE)

# Reformat protein ids
attn_aca <- attn_aca %>% mutate(V1 = gsub("_T", "_t", V1))
attn_bdi <- attn_bdi %>% mutate(V1 = sub("bdiv", "Bdiv", tolower(V1)))
ogs_long_primary_gla <- ogs_long_primary %>% filter(taxid == "5741") %>% mutate(accession_reduce = toupper(gsub("\\.", "", accession)))
attn_gla$V1 <- ogs_long_primary_gla$accession[match(attn_gla$V1, ogs_long_primary_gla$accession_reduce)]
ogs_long_primary_ath <- ogs_long_primary %>% filter(taxid == "3702") %>% mutate(accession_reduce = toupper(gsub("\\.", "", accession)))
attn_ath$V1 <- ogs_long_primary_ath$accession[match(attn_ath$V1, ogs_long_primary_ath$accession_reduce)]
ogs_long_primary_tbr <- ogs_long_primary %>% filter(taxid == "185431") %>% mutate(accession_reduce = toupper(gsub("\\.", "", accession)))
attn_tbr$V1 <- ogs_long_primary_tbr$accession[match(attn_tbr$V1, ogs_long_primary_tbr$accession_reduce)]
ogs_long_primary_lta <- ogs_long_primary %>% filter(taxid == "5689") %>% mutate(accession_reduce = toupper(gsub("\\.", "", accession)))
attn_lta$V1 <- ogs_long_primary_lta$accession[match(attn_lta$V1, ogs_long_primary_lta$accession_reduce)]

# Combine
attn_raw <- rbind(attn_human, attn_yeast, attn_aca, attn_gla, attn_ath, attn_bdi, attn_tbr, attn_lta)
colnames(attn_raw) <- c("protein_id", "AA", "alpha")

attn_raw <- attn_raw %>% group_by(protein_id) %>% mutate(position_n = row_number()) %>% mutate(position_c = max(position_n) - position_n + 1) %>% ungroup()

# Remove organelle-encoded proteins
attn <- attn_raw %>% filter(!protein_id %in% c(all_mito_accessions, all_nonmito_accessions))

## Filter to include only proteins predicted mito by DeepLoc
# Read in Retrained DeepLoc data
deeploc_results <- read.table(here("data/deeploc/predictions", "DeepLoc2.0-mito_predictions.tsv"), header=TRUE)
colnames(deeploc_results) <- c("Protein_ID", "Mitochondrion")
deeploc_results$taxid <- gsub("_.*", "", deeploc_results$Protein_ID)
deeploc_thresholds <- read.csv(here("data/deeploc", "deeploc_thresholds.csv"))
mito_threshold <- deeploc_thresholds$threshold[deeploc_thresholds$label == "Mitochondrion"]
deeploc_results$label <- "Nonmito"
deeploc_results$label[deeploc_results$Mitochondrion >= mito_threshold] <- "Mito"
deeploc_results_mito <- deeploc_results %>% filter(label == "Mito")
attn <- attn %>% filter(protein_id %in% deeploc_results_mito$Protein_ID)

# Filter to include only MitoCarta proteins
attn <- attn %>% filter(protein_id %in% gold_gene_accession_OG_id_df$gene_accession)

# Add OG IDs
attn <- merge(attn, ogs_long_primary, by.x="protein_id", by.y="accession")

# Add TargetP
attn$targetp_prediction <- targetp_results$Prediction[match(attn$protein_id, targetp_results$ID)]
attn$targetp_mTP_prob <- targetp_results$mTP[match(attn$protein_id, targetp_results$ID)]
attn$MTS_label <- "No MTS"
attn$MTS_label[attn$targetp_mTP_prob > possible_mts_threshold] <- "Possible MTS"
attn$MTS_label[attn$targetp_prediction == "mTP"] <- "MTS"
attn$MTS_label <- factor(attn$MTS_label, levels = c("No MTS", "Possible MTS", "MTS"))

## Map to human MitoCarta
attn_mitocarta <- merge(attn, human_mitocarta, by.x="protein_id", by.y="gene_accession", all.x=TRUE)
attn_mitocarta <- attn_mitocarta %>% group_by(Orthogroup) %>% mutate(MitoCarta3.0_SubMitoLocalization = paste0(unique(MitoCarta3.0_SubMitoLocalization[!is.na(MitoCarta3.0_SubMitoLocalization)]), collapse="|"), human_homologs = paste(unique(Symbol.x[!is.na(Symbol.x)]), collapse=","))

# Plot by submitochondrial localization
attn_mitocarta_summary <- attn_mitocarta %>% filter(MitoCarta3.0_SubMitoLocalization %in% selected_submito_compartments) %>% group_by(MitoCarta3.0_SubMitoLocalization, taxid) %>% mutate(alpha_normalized = alpha / sum(alpha)) %>% ungroup() %>% group_by(protein_id) %>% mutate(position_n_normalized = position_n / max(position_n)) %>% ungroup()
attn_mitocarta_summary_counts <- attn_mitocarta_summary %>% group_by(MitoCarta3.0_SubMitoLocalization, taxid) %>% summarize(n_proteins = length(unique(protein_id)))
attn_mitocarta_summary$submitolabel <- factor(attn_mitocarta_summary$MitoCarta3.0_SubMitoLocalization, levels = selected_submito_compartments)
attn_mitocarta_summary$taxid <- factor(attn_mitocarta_summary$taxid, levels = selected_taxids)
pdf('attention_mitocarta_mitotol_submito_deeploc.mito.pdf', width = 20, height = 10)
ggplot(data = attn_mitocarta_summary, aes(x = position_n_normalized, weight = alpha_normalized, fill = MTS_label)) + geom_histogram(binwidth = bin_width) + facet_wrap(~ submitolabel + taxid, scales = "fixed", ncol = 8, drop = FALSE) + scale_fill_manual(values = c("MTS" = "black", "Possible MTS" = "#999999", "No MTS" = "#00A36C"))
dev.off()

