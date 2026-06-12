library(tidyverse)
library(ape)
library(phytools)
library(TreeTools)
library(ggplot2)
library(ggrepel)
library(gplots)
library(Matrix)
library(ggpattern)

# Set ggplot theme
theme_set(theme_classic())

# Mitoproteome summary barplots
gold_gene_accession_OG_id_primary_df <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2026.04.05.tsv"), sep="\t", header=TRUE)

# Read in OGs for 673 euks
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk673spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")

mtdna_encoded_protein_count <- read.csv(here("data", "comparative", "mtdna_encoded_protein_counts.csv"))

# Get PFAM annotations
pfam_seprows <- read.table(here("data", "annotation", "pfam", "euk673spp_prokgroups_pfam_list.filtered_seprows.tsv"), sep="\t", header=TRUE)
gold_gene_accession_OG_id_df_pfam <- merge(gold_gene_accession_OG_id_primary_df, pfam_seprows, by.x = "gene_accession", by.y = "Accession", all.x = TRUE)

completed_mitoproteomes_species_list <- c("5741", "1257118", "9606", "559292", "32595", "3702", "185431", "5689")

total_protein_count <- ogs_long %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% group_by(taxid) %>% summarize(n_total_proteins = length(unique(accession)))
mito_protein_count <- gold_gene_accession_OG_id_primary_df %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% group_by(taxid) %>% summarize(n_mito_proteins = length(unique(gene_accession)))

protein_count <- merge(total_protein_count, mito_protein_count, by="taxid")
protein_count <- merge(protein_count, mtdna_encoded_protein_count, by="taxid")

mito_pfam_count <- gold_gene_accession_OG_id_df_pfam %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% filter(!is.na(PfamList)) %>% group_by(taxid) %>% summarize(n_unique_pfam = length(unique(PfamList)), n_mito_proteins_with_pfam = length(unique(gene_accession)))

protein_count <- merge(protein_count, mito_pfam_count, by="taxid")

protein_count <- protein_count %>% mutate(fraction_mito = n_mito_proteins / n_total_proteins, fraction_mito_with_pfam = n_mito_proteins_with_pfam / n_mito_proteins)

protein_count_long <- protein_count %>% pivot_longer(-taxid)
selected_vars <- c("n_total_proteins", "n_mtdna_encoded_proteins", "n_mito_proteins", "fraction_mito", "fraction_mito_with_pfam")
protein_count_long <- protein_count_long %>% filter(name %in% selected_vars)
protein_count_long$name <- factor(protein_count_long$name, levels = selected_vars)
protein_count_long$taxid <- factor(protein_count_long$taxid, levels = completed_mitoproteomes_species_list)

pdf('mitoproteome_summary_stats_barplot.pdf', width=8, height=5)
ggplot(data = protein_count_long, aes(x = value, y = taxid, label = round(value, digits=2))) + geom_col(fill = "black") + geom_text(hjust=-0.1) + facet_wrap(name ~ ., scales = "free_x", nrow = 1) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0.4)))
dev.off()



## Plot shared protein counts

# Pie chart for abSENSE across single-clade OGs
homology_power_agg <- read.table(here("data/abSENSE_HMM", "species_tree_1", "absense_results.tsv"), sep="\t", header=TRUE)
colnames(homology_power_agg) <- c("PhROG_id", "human_accessions", "L", "R", "OG_mrca_label", "OG_outgroup_mrca_label", "outgroup_kingdom_species", "fraction_of_detectable_outgroup_kingdoms", "probability_of_detection_in_any_outgroup_species", "hmm_length", "r_squared", "comments")
homology_power_agg <- homology_power_agg %>% rowwise() %>% mutate(OG_id = gsub("_Node.*", "", PhROG_id))
# Power to detect in any outgroup kingdoms
OG_ids_detectable <- unique(homology_power_agg$OG_id[which(homology_power_agg$fraction_of_detectable_outgroup_kingdoms > 0 | homology_power_agg$OG_mrca_label == "Node34_Eukaryota")])

completed_mitoproteomes_species_list <- c("5689", "185431", "3702", "32595", "9606", "559292", "1257118", "5741")

total_protein_count_shared <- ogs_long %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% group_by(Orthogroup) %>% mutate(n_clades_shared = sum("5741" %in% taxid, "1257118" %in% taxid, "9606" %in% taxid, "559292" %in% taxid, "32595" %in% taxid, "3702" %in% taxid, any(c("185431", "5689") %in% taxid))) %>% filter(BOOL_PRIMARY_OG)
mito_protein_count_shared <- ogs_long %>% filter(accession %in% gold_gene_accession_OG_id_primary_df$gene_accession) %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% group_by(Orthogroup) %>% mutate(n_clades_shared = sum("5741" %in% taxid, "1257118" %in% taxid, "9606" %in% taxid, "559292" %in% taxid, "32595" %in% taxid, "3702" %in% taxid, any(c("185431", "5689") %in% taxid))) %>% filter(BOOL_PRIMARY_OG)

total_protein_count_shared$label <- total_protein_count_shared$n_clades_shared
mito_protein_count_shared$label <- mito_protein_count_shared$n_clades_shared
mito_protein_count_shared$label[mito_protein_count_shared$n_clades_shared == 1 & mito_protein_count_shared$Orthogroup %in% OG_ids_detectable] <- "1:Powered"
mito_protein_count_shared$label[mito_protein_count_shared$n_clades_shared == 1 & !mito_protein_count_shared$Orthogroup %in% OG_ids_detectable] <- "1:Not powered"

total_protein_count_shared_summary <- total_protein_count_shared %>% group_by(taxid, label) %>% summarize(value = length(unique(accession)))
total_protein_fraction_shared_summary <- total_protein_count_shared %>% group_by(taxid) %>% mutate(n_total_proteins = length(unique(accession))) %>% group_by(taxid, label, n_total_proteins) %>% summarize(n_proteins = length(unique(accession))) %>% mutate(value = n_proteins / n_total_proteins)
mito_protein_count_shared_summary <- mito_protein_count_shared %>% group_by(taxid, label) %>% summarize(value = length(unique(accession)))
mito_protein_fraction_shared_summary <- mito_protein_count_shared %>% group_by(taxid) %>% mutate(n_total_proteins = length(unique(accession))) %>% group_by(taxid, label, n_total_proteins) %>% summarize(n_proteins = length(unique(accession))) %>% mutate(value = n_proteins / n_total_proteins)

total_protein_count_shared_summary$name <- "Total proteins"
total_protein_fraction_shared_summary$name <- "Fraction total proteins"
mito_protein_count_shared_summary$name <- "Mito proteins"
mito_protein_fraction_shared_summary$name <- "Fraction mito proteins"

selected_columns <- c("name", "taxid", "label", "value")
protein_count_summary <- rbind(mito_protein_count_shared_summary[,selected_columns], mito_protein_fraction_shared_summary[,selected_columns])

protein_count_summary$name <- factor(protein_count_summary$name, levels = c("Mito proteins", "Fraction mito proteins"))
protein_count_summary$taxid <- factor(protein_count_summary$taxid, levels = rev(completed_mitoproteomes_species_list))
protein_count_summary$label <- as.factor(protein_count_summary$label)
protein_count_summary$species_name <- uniprot_proteomes_tax$ScientificName[match(protein_count_summary$taxid, uniprot_proteomes_tax$TaxId)]
protein_count_summary$species_name <- factor(protein_count_summary$species_name, levels = uniprot_proteomes_tax$ScientificName[match(rev(completed_mitoproteomes_species_list), uniprot_proteomes_tax$TaxId)])

library(colorspace)
base_color <- "#0047AB"
# Generate lighter shades
lighter_1 <- lighten(base_color, 0.8)
lighter_2 <- lighten(base_color, 0.6)
lighter_3 <- lighten(base_color, 0.4)
lighter_4 <- lighten(base_color, 0.2)
lighter_5 <- lighten(base_color, 0.1)
colors <- c("grey90", "grey90", lighter_1, lighter_2, lighter_3, lighter_4, lighter_5, base_color)

pdf('mitoproteome_summary_stats_barplot_proteins_shared.pdf', width=6, height=4)
ggplot(data = protein_count_summary, aes(x = value, y = species_name, group = label, fill = label, pattern = label)) + geom_col_pattern(width = 0.75, linewidth = 0.4, color = "black", pattern_colour = "white", pattern_fill = "white", pattern_size = 0.1, pattern_angle = 45) + scale_fill_manual(values = colors) + scale_pattern_manual(values = c("stripe", "none", "none", "none", "none", "none", "none", "none")) + facet_wrap(name ~ ., scales = "free_x", nrow = 1)
dev.off()


# With pattern
mito_protein_count_shared_pie_summary <- ogs_long %>% filter(accession %in% gold_gene_accession_OG_id_df$gene_accession) %>% filter(taxid %in% completed_mitoproteomes_species_list) %>% group_by(Orthogroup) %>% summarize(n_clades_shared = sum("5741" %in% taxid, "1257118" %in% taxid, "9606" %in% taxid, "559292" %in% taxid, "32595" %in% taxid, "3702" %in% taxid, any(c("185431", "5689") %in% taxid)))
mito_protein_count_shared_pie_summary$label <- mito_protein_count_shared_pie_summary$n_clades_shared
mito_protein_count_shared_pie_summary$label[mito_protein_count_shared_pie_summary$n_clades_shared == 1 & !mito_protein_count_shared_pie_summary$Orthogroup %in% OG_ids_detectable] <- "1:Not powered"
mito_protein_count_shared_pie_summary$label[mito_protein_count_shared_pie_summary$n_clades_shared == 1 & mito_protein_count_shared_pie_summary$Orthogroup %in% OG_ids_detectable] <- "1:Powered"
mito_protein_count_shared_pie_summary_counts <- mito_protein_count_shared_pie_summary %>% group_by(label) %>% summarize(n = n())
pdf('mitoproteome_summary_stats_pie_proteins_shared_absense.pdf', width=4, height=4)
ggplot(mito_protein_count_shared_pie_summary_counts, aes(x = "", y = n, fill = label, pattern = label)) +
  geom_bar_pattern(stat = "identity", width = 1, linewidth = 0.4, color = "black", pattern_colour = "white", pattern_fill = "white", pattern_size = 0.1, pattern_angle = 45) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colors) +
  scale_pattern_manual(values = c("stripe", "none", "none", "none", "none", "none", "none", "none")) +
  theme_void()
dev.off()




## Compare mitoproteomes with Venn and Upset

# Read in OGs for 203 euks
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")

uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

selected_taxids <- c("5689", "185431", "3702", "32595", "9606", "559292", "1257118", "5741")

# Include non-primary OGs
gold_gene_accession_OG_id_df <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_2026.04.05.tsv"), sep="\t", header=TRUE)
gold_gene_accession_OG_id_df$taxid <- as.character(gold_gene_accession_OG_id_df$taxid)
gold_gene_accession_OG_id_subset_df <- gold_gene_accession_OG_id_df[gold_gene_accession_OG_id_df$taxid %in% selected_taxids,]
gold_gene_accession_OG_id_subset_df %>% filter(OG_id != "") %>% group_by(taxid) %>% summarize(n = length(unique(gene_accession)))
gold_gene_accession_OG_id_subset_df <- gold_gene_accession_OG_id_subset_df %>% filter(OG_id != "")

gold_gene_accession_OG_id_subset_df_fusion_counts <- gold_gene_accession_OG_id_subset_df %>% mutate(multiple_OGs = grepl(",", OG_id, fixed=TRUE)) %>% group_by(taxid) %>% summarize(n_mito_proteins = n(), n_mito_proteins_in_multiple_OGs = sum(multiple_OGs)) %>% mutate(n_mito_proteins_in_single_OG = n_mito_proteins - n_mito_proteins_in_multiple_OGs, fraction_mito_proteins_in_multiple_OGs = round(n_mito_proteins_in_multiple_OGs / n_mito_proteins, digits=3))
gold_gene_accession_OG_id_subset_df_fusion_counts$species <- uniprot_proteomes_tax$ScientificName[match(gold_gene_accession_OG_id_subset_df_fusion_counts$taxid, uniprot_proteomes_tax$TaxId)]

gold_gene_accession_OG_id_subset_df <- gold_gene_accession_OG_id_subset_df %>% separate_rows(OG_id, sep=",")
gold_gene_accession_OG_id_subset_human <- gold_gene_accession_OG_id_subset_df %>% filter(taxid == "9606")

# Get primary OGs
gold_gene_accession_OG_id_primary_df$taxid <- as.character(gold_gene_accession_OG_id_primary_df$taxid)
gold_gene_accession_OG_id_subset_primary_df <- gold_gene_accession_OG_id_primary_df[gold_gene_accession_OG_id_primary_df$taxid %in% selected_taxids,]
gold_gene_accession_OG_id_subset_primary_df %>% filter(OG_id != "") %>% group_by(taxid) %>% summarize(n = length(unique(gene_accession)))
gold_gene_accession_OG_id_subset_primary_df <- gold_gene_accession_OG_id_subset_primary_df %>% filter(OG_id != "")
gold_gene_accession_OG_id_subset_primary_human <- gold_gene_accession_OG_id_subset_primary_df %>% filter(taxid == "9606")

gold_mito_OG_counts_per_species <- gold_gene_accession_OG_id_subset_df %>% group_by(taxid) %>% summarize(n_OGs = length(unique(OG_id)))

human_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "9606"]
yeast_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "559292"]
tbr_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "185431"]
lta_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "5689"]
aca_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "1257118"]
ath_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "3702"]
gla_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "5741"]
bdi_mito_OGs <- gold_gene_accession_OG_id_subset_df$OG_id[gold_gene_accession_OG_id_subset_df$taxid == "32595"]

mito_sets <- list(`Human` = unlist(strsplit(human_mito_OGs, split=",")), `Yeast` = unlist(strsplit(yeast_mito_OGs, split=",")), `Acanthamoeba` = unlist(strsplit(aca_mito_OGs, split=",")), `Giardia` = unlist(strsplit(gla_mito_OGs, split=",")), `Trypanosoma` = unlist(strsplit(tbr_mito_OGs, split=",")), `Leishmania` = unlist(strsplit(lta_mito_OGs, split=",")), `Arabidopsis` = unlist(strsplit(ath_mito_OGs, split=",")), `Babesia` = unlist(strsplit(bdi_mito_OGs, split=",")))
mito_OG_ids <- unique(unlist(mito_sets))


## UpSet with absense annotated stacked barplot
# Get absense predictions
mito_OG_ids_detectable <- mito_OG_ids[mito_OG_ids %in% OG_ids_detectable]

# Get absense predictions for human
homology_power_per_species <- read.table(here("data/abSENSE_HMM", "species_tree_1", "absense_results_per_species.tsv"), sep="\t", header=FALSE)
colnames(homology_power_per_species) <- c("PhROG_id", "taxid", "probability_of_detection")
homology_power_per_species <- homology_power_per_species %>% rowwise() %>% mutate(OG_id = gsub("_Node.*", "", PhROG_id))
homology_power_human <- homology_power_per_species %>% filter(taxid == "9606")
# Mark homologs as probability of detection = 1
ogs_long_select <- ogs_long %>% filter(taxid %in% "9606") %>% filter(Orthogroup %in% mito_OG_ids)
homology_power_human <- homology_power_human %>% mutate(bool_has_homolog = taxid %in% ogs_long_select$taxid[ogs_long_select$Orthogroup == OG_id])
homology_power_human_detectable <- homology_power_human %>% filter(probability_of_detection > 0.95 | bool_has_homolog)
mito_OG_ids_detectable_human <- unique(homology_power_human_detectable$OG_id)

gold_gene_accession_OG_id_subset_df_seprows <- gold_gene_accession_OG_id_subset_df %>% separate_rows(OG_id, sep=",")
mito_OG_ids <- unique(gold_gene_accession_OG_id_subset_df_seprows$OG_id)
mito_OGs_binary_mat <- matrix(data = 0, nrow = length(mito_OG_ids), ncol = length(selected_taxids))
rownames(mito_OGs_binary_mat) <- mito_OG_ids
colnames(mito_OGs_binary_mat) <- selected_taxids
for (i in 1:length(mito_OG_ids)) {
  curr_OG_id <- mito_OG_ids[i]
  gold_gene_accession_OG_id_subset_df_curr <- gold_gene_accession_OG_id_subset_df_seprows %>% filter(OG_id == curr_OG_id)
  mito_OGs_binary_mat[i, unique(gold_gene_accession_OG_id_subset_df_curr$taxid)] <- 1
}
species_abbrev <- gsub(" .*", "", uniprot_proteomes_tax$ScientificName[match(selected_taxids, uniprot_proteomes_tax$TaxId)])
colnames(mito_OGs_binary_mat) <- species_abbrev
mito_OGs_binary <- as.data.frame(mito_OGs_binary_mat)
mito_OGs_binary$Orthogroup <- rownames(mito_OGs_binary_mat)
mito_OGs_binary$homology_detection_power <- "Not powered"
mito_OGs_binary$homology_detection_power[which(rownames(mito_OGs_binary) %in% mito_OG_ids_detectable)] <- "Powered" # completed mitoproteomes
mito_OGs_binary$homology_detection_power <- factor(mito_OGs_binary$homology_detection_power, levels = c("Not powered", "Powered"))

mito_OGs_binary$n_clades_present <- rowSums(mito_OGs_binary[,c("Homo", "Saccharomyces", "Acanthamoeba", "Giardia", "Arabidopsis", "Babesia")]) + as.numeric(rowSums(mito_OGs_binary[,c("Trypanosoma", "Leishmania")]) > 0)
mito_OGs_binary$label <- mito_OGs_binary$n_clades_present
mito_OGs_binary$label[mito_OGs_binary$n_clades_present == 1 & mito_OGs_binary$homology_detection_power == "Powered"] <- "1:Powered"
mito_OGs_binary$label[mito_OGs_binary$n_clades_present == 1 & mito_OGs_binary$homology_detection_power == "Not powered"] <- "1:Not powered"

library(colorspace)
base_color <- "#0047AB"
# Generate lighter shades
lighter_1 <- lighten(base_color, 0.8)
lighter_2 <- lighten(base_color, 0.6)
lighter_3 <- lighten(base_color, 0.4)
lighter_4 <- lighten(base_color, 0.2)
lighter_5 <- lighten(base_color, 0.1)
colors <- c("grey90", "grey90", lighter_1, lighter_2, lighter_3, lighter_4, lighter_5, base_color)

library(ComplexUpset)
pdf("mito_OGs_complexupset_detectable.euks.pdf", width = 16, height = 6)
ComplexUpset::upset(
  mito_OGs_binary,
  rev(species_abbrev),
  base_annotations=list(),
  annotations = list(
    'Intersection' = (
      ggplot(mapping = aes(x = intersection, fill = label, pattern = label))
      + geom_bar_pattern(linewidth = 0.4, color = "black", pattern_colour = "white", pattern_fill = "white", pattern_size = 0.1, pattern_angle = 45)
      + geom_text(aes(label=..count..), stat='count', vjust=-0.5, size=3)
      + scale_fill_manual(values = colors)
      + scale_pattern_manual(
        values = c("stripe", "none", "none", "none", "none", "none", "none", "none"),
        guide = 'none'
      )
    ) + ylim(0,1000) + ylab("# shared orthogroups") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  ),
  set_sizes=upset_set_size(
    geom=geom_bar(
      stat='count',
      fill='black'
    )
  ) + geom_text(aes(label=..count..), hjust=1.1, size=3, stat='count') + ylim(1500,0) + ylab("# mitochondrial orthogroups") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()),
  min_size=10, # set minimum size
  sort_sets=FALSE,
  sort_intersections_by=c('degree', 'cardinality'),
  guides='over' # moves legends over the set sizes
)
dev.off()

# Map to mitopathways for manual labeling
human_only_OG_ids <- mito_sets$Human
human_only_OG_ids <- unique(human_only_OG_ids[!human_only_OG_ids %in% Reduce(union, mito_sets[c("Yeast", "Acanthamoeba", "Giardia", "Leishmania", "Trypanosoma", "Arabidopsis", "Babesia")])])
gold_gene_accession_OG_id_subset_select <- gold_gene_accession_OG_id_subset_df %>% filter(OG_id %in% human_only_OG_ids)

mitopathways_human_id2training_mito <- read.table(here("data", "annotation", "human_mitocarta", "human_mitopathways_uniprot.tsv"), sep="\t", header=TRUE)
mitopathways_human_id2training_mito$accession <- paste0("9606_", mitopathways_human_id2training_mito$Entry)
gold_gene_accession_OG_id_df_subset_human <- merge(gold_gene_accession_OG_id_subset_select, mitopathways_human_id2training_mito, by.x="gene_accession", by.y="accession")
mitopathways_human_id2training_mito_counts <- mitopathways_human_id2training_mito %>% group_by(MitoPathway) %>% summarize(n_total_members = length(unique(Gene)))
gold_gene_accession_OG_id_df_subset_counts <- gold_gene_accession_OG_id_df_subset_human %>% group_by(MitoPathway) %>% summarize(n_proteins = length(unique(gene_accession)))
gold_gene_accession_OG_id_df_subset_mitopathways <- merge(gold_gene_accession_OG_id_df_subset_counts, mitopathways_human_id2training_mito_counts, by="MitoPathway")
gold_gene_accession_OG_id_df_subset_mitopathways <- gold_gene_accession_OG_id_df_subset_mitopathways %>% mutate(fraction_mitopathway = n_proteins / n_total_members)
gold_gene_accession_OG_id_df_subset_mitopathways <- gold_gene_accession_OG_id_df_subset_mitopathways
gold_gene_accession_OG_id_df_subset_mitopathways <- gold_gene_accession_OG_id_df_subset_mitopathways[order(gold_gene_accession_OG_id_df_subset_mitopathways$fraction_mitopathway, decreasing=TRUE),]

## Pie charts
# Map PFAM annotations
gold_gene_accession_OG_id_subset_select <- gold_gene_accession_OG_id_subset_df %>% filter(taxid %in% selected_taxids) %>% separate_rows(OG_id, sep=",")
gold_gene_accession_OG_id_df_pfam <- merge(gold_gene_accession_OG_id_subset_select, pfam_seprows, by.x = "gene_accession", by.y = "Accession", all.x = TRUE)

# Add homology power
gold_gene_accession_OG_id_df_pfam$absense <- "Not powered"
gold_gene_accession_OG_id_df_pfam$absense[gold_gene_accession_OG_id_df_pfam$OG_id %in% mito_OG_ids_detectable] <- "Powered"
gold_gene_accession_OG_id_df_pfam_summary <- gold_gene_accession_OG_id_df_pfam %>% group_by(OG_id) %>% mutate(n_proteins = length(unique(gene_accession))) %>% group_by(OG_id, PfamList) %>% summarize(absense = unique(absense), n_proteins = unique(n_proteins), n_proteins_with_pfam = n()) %>% mutate(fraction_proteins_with_pfam = n_proteins_with_pfam / n_proteins)
gold_gene_accession_OG_id_df_pfam_summary <- gold_gene_accession_OG_id_df_pfam_summary %>% group_by(OG_id) %>% summarize(absense = unique(absense), majority_pfam = paste0(sort(PfamList[fraction_proteins_with_pfam >= 0.5]), collapse=","), any_pfam = paste0(sort(PfamList), collapse=","))
gold_gene_accession_OG_id_df_nopfam <- gold_gene_accession_OG_id_df_pfam_summary %>% filter(any_pfam == "")

nopfam_absense_counts <- table(gold_gene_accession_OG_id_df_nopfam$absense)
pdf("pie_mito_OGs_nopfam_absense.pdf", width = 6, height = 6)
pie(nopfam_absense_counts, labels = paste0(names(nopfam_absense_counts), " (", nopfam_absense_counts, ")"))
dev.off()

gold_gene_accession_OG_id_df_majoritypfam <- gold_gene_accession_OG_id_df_pfam_summary %>% filter(majority_pfam != "")
gold_gene_accession_OG_id_df_anypfam <- gold_gene_accession_OG_id_df_pfam_summary %>% filter(any_pfam != "" & !OG_id %in% gold_gene_accession_OG_id_df_majoritypfam$OG_id)
gold_gene_accession_OG_id_df_nopfam <- gold_gene_accession_OG_id_df_pfam_summary %>% filter(any_pfam == "")
pfam_counts <- c(nrow(gold_gene_accession_OG_id_df_majoritypfam), nrow(gold_gene_accession_OG_id_df_anypfam), nrow(gold_gene_accession_OG_id_df_nopfam))
pdf("pie_mito_OGs_pfam.pdf", width = 6, height = 6)
pie(pfam_counts, labels = c(paste0("Majority Pfam", " (", nrow(gold_gene_accession_OG_id_df_majoritypfam), ")"), paste0("Any Pfam", " (",  nrow(gold_gene_accession_OG_id_df_anypfam), ")"), paste0("No Pfam", " (", nrow(gold_gene_accession_OG_id_df_nopfam), ")")))
dev.off()

gold_gene_accession_OG_id_df_pfam_summary_pie <- gold_gene_accession_OG_id_df_pfam_summary
gold_gene_accession_OG_id_df_pfam_summary_pie$label <- ""
gold_gene_accession_OG_id_df_pfam_summary_pie$label[gold_gene_accession_OG_id_df_pfam_summary_pie$absense == "Powered"] <- "No Pfam:Powered"
gold_gene_accession_OG_id_df_pfam_summary_pie$label[gold_gene_accession_OG_id_df_pfam_summary_pie$absense == "Not powered"] <- "No Pfam:Not powered"
gold_gene_accession_OG_id_df_pfam_summary_pie$label[gold_gene_accession_OG_id_df_pfam_summary_pie$majority_pfam != ""] <- "Majority Pfam"
gold_gene_accession_OG_id_df_pfam_summary_pie$label[gold_gene_accession_OG_id_df_pfam_summary_pie$any_pfam != "" & gold_gene_accession_OG_id_df_pfam_summary_pie$majority_pfam == ""] <- "Any Pfam"
gold_gene_accession_OG_id_df_pfam_summary_pie <- gold_gene_accession_OG_id_df_pfam_summary_pie
gold_gene_accession_OG_id_df_pfam_summary_pie_counts <- gold_gene_accession_OG_id_df_pfam_summary_pie %>% group_by(label) %>% summarize(n = n())
gold_gene_accession_OG_id_df_pfam_summary_pie_counts$label <- factor(gold_gene_accession_OG_id_df_pfam_summary_pie_counts$label, levels = c("Majority Pfam", "Any Pfam", "No Pfam:Powered", "No Pfam:Not powered"))
pdf("pie_mito_OGs_pfam_absense.pdf", width = 6, height = 6)
ggplot(gold_gene_accession_OG_id_df_pfam_summary_pie_counts, aes(x = "", y = n, fill = label, pattern = label)) +
  geom_bar_pattern(stat = "identity", width = 1, linewidth = 0.4, color = "black", pattern_colour = "white", pattern_fill = "white", pattern_size = 0.1, pattern_angle = 45) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#662d91", "#d58aff", "#d8d8d8", "#d8d8d8")) +
  scale_pattern_manual(values = c("none", "none", "none", "stripe")) +
  theme_void()
dev.off()

gold_gene_accession_OG_id_df_pfam_summary_human <- gold_gene_accession_OG_id_df_pfam_summary %>% filter(OG_id %in% mito_sets$Human)
gold_gene_accession_OG_id_df_majoritypfam <- gold_gene_accession_OG_id_df_pfam_summary_human %>% filter(majority_pfam != "")
gold_gene_accession_OG_id_df_anypfam <- gold_gene_accession_OG_id_df_pfam_summary_human %>% filter(any_pfam != "" & !OG_id %in% gold_gene_accession_OG_id_df_majoritypfam$OG_id)
gold_gene_accession_OG_id_df_nopfam <- gold_gene_accession_OG_id_df_pfam_summary_human %>% filter(any_pfam == "")
pfam_counts <- c(nrow(gold_gene_accession_OG_id_df_majoritypfam), nrow(gold_gene_accession_OG_id_df_anypfam), nrow(gold_gene_accession_OG_id_df_nopfam))
pdf("pie_human_mito_OGs_pfam_2026.05.11.pdf", width = 6, height = 6)
pie(pfam_counts, labels = c(paste0("Majority Pfam", " (", nrow(gold_gene_accession_OG_id_df_majoritypfam), ")"), paste0("Any Pfam", " (",  nrow(gold_gene_accession_OG_id_df_anypfam), ")"), paste0("No Pfam", " (", nrow(gold_gene_accession_OG_id_df_nopfam), ")")))
dev.off()

gold_gene_accession_OG_id_subset_counts <- gold_gene_accession_OG_id_subset_df %>% group_by(OG_id) %>% summarize(n_species = length(unique(taxid)), n_proteins = length(unique(gene_accession)))
pdf("barplot_species_per_MitoCarta_OG.pdf", width = 6, height = 4)
ggplot(data = gold_gene_accession_OG_id_subset_counts, aes(x = n_species)) + geom_bar(fill="black") + xlab("Species with MitoCarta homolog") + ylab("MitoCarta Orthogroups")
dev.off()
pdf("barplot_proteins_per_MitoCarta_OG.pdf", width = 6, height = 4)
ggplot(data = gold_gene_accession_OG_id_subset_counts, aes(x = n_proteins)) + geom_bar() + xlab("MitoCarta proteins") + ylab("MitoCarta Orthogroups")
dev.off()


gold_gene_accession_OG_id_subset_primary_human_yeast <- gold_gene_accession_OG_id_subset_human %>% filter(OG_id %in% mito_sets$Yeast)
yeast_OG_ids <- unique(gold_gene_accession_OG_id_subset_primary_human_yeast$OG_id)
other_OG_ids <- Reduce(union, mito_sets[c("Acanthamoeba", "Giardia", "Leishmania", "Trypanosoma", "Arabidopsis", "Babesia")])
gold_gene_accession_OG_id_subset_primary_human_anyother <- gold_gene_accession_OG_id_subset_human %>% filter(OG_id %in% other_OG_ids)
shared_other_OG_ids <- unique(gold_gene_accession_OG_id_subset_primary_human_anyother$OG_id)
yeast_and_other <- intersect(yeast_OG_ids, shared_other_OG_ids)
yeast_only <- yeast_OG_ids[!yeast_OG_ids %in% other_OG_ids]
other_only <- shared_other_OG_ids[!shared_other_OG_ids %in% yeast_OG_ids]
not_shared <- unique(gold_gene_accession_OG_id_subset_human$OG_id[!gold_gene_accession_OG_id_subset_human$OG_id %in% c(yeast_OG_ids, shared_other_OG_ids)])
OG_counts <- c(length(yeast_and_other), length(other_only), length(yeast_only), length(not_shared))
pdf("pie_human_mito_OGs_shared.pdf", width = 6, height = 6)
pie(OG_counts, labels = c(paste0("Yeast and other species", " (", OG_counts[1], ")"), paste0("Other species", " (",  OG_counts[2], ")"), paste0("Yeast", " (", OG_counts[3], ")"), paste0("Not shared", " (", OG_counts[4], ")")))
dev.off()


# Map to human homologs and mitopathways
human_annot <- read.table(here("data", "orthogroups", "idmapping", "hsa_geneid_to_symbol.tsv"), sep="\t")
human_annot$accession <- paste0("9606_", human_annot$V1)
human_mito_uniprot_ids <- read.table(here("data", "annotation", "human_mitocarta", "human.mito.id_uniprot.txt"), header=TRUE)
human_id2training <- read.table(here("data", "annotation", "human_mitocarta", "human.id2training_uniprot_symbols.txt"), sep="\t", header=TRUE)
human_id2training_mito <- human_id2training
human_id2training_mito <- human_id2training_mito[order(human_id2training_mito$Symbol),]
## Map to human MitoPathways
mitopathways <- read.table(here("data", "annotation", "human_mitocarta", "human.path2gene.txt"), sep="\t", header=TRUE)
mitopathways_human_id2training_mito <- merge(mitopathways, human_id2training_mito, by.x="Gene", by.y="Symbol")
mitopathways_human_id2training_mito$gene_accession <- paste0("9606_", mitopathways_human_id2training_mito$Entry)
gold_gene_accession_OG_id_subset_primary_mitopathway <- merge(gold_gene_accession_OG_id_subset_human, mitopathways_human_id2training_mito, by="gene_accession")

gold_gene_accession_OG_id_subset_primary_human_yeast <- gold_gene_accession_OG_id_subset_human %>% filter(OG_id %in% mito_sets$Yeast)
yeast_OG_ids <- unique(gold_gene_accession_OG_id_subset_primary_human_yeast$OG_id)
other_OG_ids <- Reduce(union, mito_sets[c("Acanthamoeba", "Giardia", "Leishmania", "Trypanosoma", "Arabidopsis", "Babesia")])
gold_gene_accession_OG_id_subset_primary_human_anyother <- gold_gene_accession_OG_id_subset_human %>% filter(OG_id %in% other_OG_ids)
shared_other_OG_ids <- unique(gold_gene_accession_OG_id_subset_primary_human_anyother$OG_id)
yeast_and_other <- intersect(yeast_OG_ids, shared_other_OG_ids)
yeast_only <- yeast_OG_ids[!yeast_OG_ids %in% other_OG_ids]
other_only <- shared_other_OG_ids[!shared_other_OG_ids %in% yeast_OG_ids]
not_shared <- unique(gold_gene_accession_OG_id_subset_human$OG_id[!gold_gene_accession_OG_id_subset_human$OG_id %in% c(yeast_OG_ids, shared_other_OG_ids)])
OG_counts <- c(length(yeast_and_other), length(other_only), length(yeast_only), length(not_shared))
gold_gene_accession_OG_id_subset_primary_mitopathway$label <- "Neither"
gold_gene_accession_OG_id_subset_primary_mitopathway$label[gold_gene_accession_OG_id_subset_primary_mitopathway$OG_id %in% yeast_and_other] <- "Yeast and other species"
gold_gene_accession_OG_id_subset_primary_mitopathway$label[gold_gene_accession_OG_id_subset_primary_mitopathway$OG_id %in% yeast_only] <- "Yeast"
gold_gene_accession_OG_id_subset_primary_mitopathway$label[gold_gene_accession_OG_id_subset_primary_mitopathway$OG_id %in% other_only] <- "Other species"
gold_gene_accession_OG_id_subset_mitopathway_counts <- gold_gene_accession_OG_id_subset_primary_mitopathway %>% group_by(MitoPathway) %>% mutate(n_orthogroups_total = length(unique(OG_id))) %>% group_by(MitoPathway, label, n_orthogroups_total) %>% summarize(n_orthogroups = length(unique(OG_id))) %>% mutate(fraction = n_orthogroups / n_orthogroups_total)

selected_mitopathways <- c("Apoptosis", "Mitophagy", "Fission", "Calcium_uniporter", "Small_molecule_transport", "CI_subunits", "CII_subunits", "CIII_subunits", "CIV_subunits", "CV_subunits", "Mitochondrial_ribosome", "mtRNA_metabolism", "mtDNA_maintenance", "Fe−S_cluster_biosynthesis", "Protein_import_and_sorting")
gold_gene_accession_OG_id_subset_mitopathway_counts <- gold_gene_accession_OG_id_subset_mitopathway_counts %>% filter(MitoPathway %in% selected_mitopathways)
gold_gene_accession_OG_id_subset_mitopathway_counts <- gold_gene_accession_OG_id_subset_mitopathway_counts %>% mutate(MitoPathway = paste0(MitoPathway, " (", n_orthogroups_total, ")"))
gold_gene_accession_OG_id_subset_mitopathway_counts_for_ordering <- gold_gene_accession_OG_id_subset_mitopathway_counts %>% filter(label != "Neither")
gold_gene_accession_OG_id_subset_mitopathway_counts_for_ordering <- gold_gene_accession_OG_id_subset_mitopathway_counts_for_ordering %>% group_by(MitoPathway) %>% summarize(fraction = sum(fraction))
mitopathway_order <- gold_gene_accession_OG_id_subset_mitopathway_counts_for_ordering$MitoPathway[order(gold_gene_accession_OG_id_subset_mitopathway_counts_for_ordering$fraction, decreasing=FALSE)]
gold_gene_accession_OG_id_subset_mitopathway_counts$MitoPathway <- factor(gold_gene_accession_OG_id_subset_mitopathway_counts$MitoPathway, levels = mitopathway_order)
gold_gene_accession_OG_id_subset_mitopathway_counts$label <- factor(gold_gene_accession_OG_id_subset_mitopathway_counts$label, levels = c("Neither", "Yeast", "Other species", "Yeast and other species"))

pdf("column_human_mitopathways_OGs_shared.pdf", width = 4, height = 8)
ggplot(data = gold_gene_accession_OG_id_subset_mitopathway_counts, aes(x = MitoPathway, y = fraction, fill = label)) + geom_col() + scale_fill_manual(values = c("#e6e6e6", "#b3b3b3", "#808080", "black")) + coord_flip() + theme(legend.position = "bottom")
dev.off()


## Pan-pathogen targets using PhROGs
selected_taxids <- c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595") # completed mitoproteomes

species_tree <- read.tree(here("data/species_phylogeny/processed_species_tree", paste0("species_tree_1.nwk")))

phrogs_long_eukaryota_parent <- read.table(here("data/phylogenetically_resolved_orthogroups", "species_tree_1", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long_eukaryota_parent) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_parent <- phrogs_long_eukaryota_parent %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% filter(!duplicated(protein_id))

gold_gene_accession_OG_id_df$taxid <- as.character(gold_gene_accession_OG_id_df$taxid)
gold_gene_accession_OG_id_subset_df <- gold_gene_accession_OG_id_df %>% filter(taxid %in% selected_taxids)

phrogs_long_eukaryota_parent_gold <- phrogs_long_eukaryota_parent %>% filter(protein_id %in% gold_gene_accession_OG_id_subset_df$gene_accession)

human_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "9606"]
yeast_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "559292"]
tbr_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "185431"]
lta_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "5689"]
aca_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "1257118"]
ath_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "3702"]
gla_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "5741"]
bdi_mito_OGs <- phrogs_long_eukaryota_parent_gold$PROG_id[phrogs_long_eukaryota_parent_gold$taxid == "32595"]

mito_sets <- list(`Human` = unlist(strsplit(human_mito_OGs, split=",")), `Yeast` = unlist(strsplit(yeast_mito_OGs, split=",")), `Acanthamoeba` = unlist(strsplit(aca_mito_OGs, split=",")), `Giardia` = unlist(strsplit(gla_mito_OGs, split=",")), `Trypanosoma` = unlist(strsplit(tbr_mito_OGs, split=",")), `Leishmania` = unlist(strsplit(lta_mito_OGs, split=",")), `Arabidopsis` = unlist(strsplit(ath_mito_OGs, split=",")), `Babesia` = unlist(strsplit(bdi_mito_OGs, split=",")))
mito_PhROG_ids <- unique(unlist(mito_sets))
length(mito_PhROG_ids)

# Get absense predictions for human
homology_power_human <- homology_power_per_species %>% filter(taxid == "9606")
# Mark homologs as probability of detection = 1
phrogs_long_eukaryota_parent_select <- phrogs_long_eukaryota_parent %>% filter(taxid %in% "9606") %>% filter(PROG_id %in% mito_PhROG_ids)
homology_power_human <- homology_power_human %>% rowwise() %>% mutate(bool_has_homolog = taxid %in% phrogs_long_eukaryota_parent_select$taxid[phrogs_long_eukaryota_parent_select$PROG_id == PhROG_id])
homology_power_human_detectable <- homology_power_human %>% filter(probability_of_detection > 0.95 | bool_has_homolog)
mito_PhROG_ids_detectable_human <- unique(homology_power_human_detectable$PhROG_id)

mito_OGs_binary_mat <- matrix(data = 0, nrow = length(mito_PhROG_ids), ncol = length(selected_taxids))
rownames(mito_OGs_binary_mat) <- mito_PhROG_ids
colnames(mito_OGs_binary_mat) <- selected_taxids
for (i in 1:length(mito_PhROG_ids)) {
  curr_OG_id <- mito_PhROG_ids[i]
  phrogs_long_eukaryota_parent_gold_curr <- phrogs_long_eukaryota_parent_gold %>% filter(PROG_id == curr_OG_id)
  mito_OGs_binary_mat[i, unique(phrogs_long_eukaryota_parent_gold_curr$taxid)] <- 1
}
species_abbrev <- gsub(" .*", "", uniprot_proteomes_tax$ScientificName[match(selected_taxids, uniprot_proteomes_tax$TaxId)])
colnames(mito_OGs_binary_mat) <- species_abbrev
mito_OGs_binary <- as.data.frame(mito_OGs_binary_mat)
mito_OGs_binary$Orthogroup <- rownames(mito_OGs_binary_mat)
mito_OGs_binary$homology_detection_power <- "Not powered"
mito_OGs_binary$homology_detection_power[which(rownames(mito_OGs_binary) %in% mito_PhROG_ids_detectable_human)] <- "Powered" # parasite drug targets


## Find parasite targets
selected_taxids <- c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595") # completed mitoproteomes

# Determine what genes are present in intersections
kinetoplast_OG_ids <- intersect(mito_sets$Leishmania, mito_sets$Trypanosoma)
parasite_OG_ids <- intersect(mito_sets$Acanthamoeba, kinetoplast_OG_ids)
parasite_OG_ids <- unique(c(kinetoplast_OG_ids, mito_sets$Acanthamoeba, mito_sets$Babesia))
length(parasite_OG_ids) # number of parasite mitochondrial PhROGs
parasite_OG_ids <- parasite_OG_ids[!parasite_OG_ids %in% c(mito_sets$Human)]
gold_gene_accession_OG_id_df_parasite <- phrogs_long_eukaryota_parent_gold[phrogs_long_eukaryota_parent_gold$PROG_id %in% parasite_OG_ids,]

gold_gene_accession_OG_id_df_parasite$gene_accession <- gold_gene_accession_OG_id_df_parasite$protein_id

# Get presence/absence in each kingdom
gold_gene_accession_OG_id_df_parasite$present_Aca_mito <- 0
gold_gene_accession_OG_id_df_parasite$present_Aca_mito[gold_gene_accession_OG_id_df_parasite$PROG_id %in% mito_sets$Acanthamoeba] <- 1
gold_gene_accession_OG_id_df_parasite$present_kinetoplastid_mito <- 0
gold_gene_accession_OG_id_df_parasite$present_kinetoplastid_mito[gold_gene_accession_OG_id_df_parasite$PROG_id %in% kinetoplast_OG_ids] <- 1
gold_gene_accession_OG_id_df_parasite$present_Bdi_mito <- 0
gold_gene_accession_OG_id_df_parasite$present_Bdi_mito[gold_gene_accession_OG_id_df_parasite$PROG_id %in% mito_sets$Babesia] <- 1
gold_gene_accession_OG_id_df_parasite$present_n_parasite_clades <- gold_gene_accession_OG_id_df_parasite$present_Aca_mito + gold_gene_accession_OG_id_df_parasite$present_kinetoplastid_mito + gold_gene_accession_OG_id_df_parasite$present_Bdi_mito

# Mark whether any human homolog exists
phrogs_long_eukaryota_parent_human_collapse <- phrogs_long_eukaryota_parent %>% filter(taxid == "9606") %>% group_by(PROG_id) %>% summarize(human_protein_ids = paste0(protein_id, collapse=", "))
gold_gene_accession_OG_id_df_parasite$human_nonmito_homologs <- phrogs_long_eukaryota_parent_human_collapse$human_protein_ids[match(gold_gene_accession_OG_id_df_parasite$PROG_id, phrogs_long_eukaryota_parent_human_collapse$PROG_id)]
gold_gene_accession_OG_id_df_parasite$human_nonmito_homologs[is.na(gold_gene_accession_OG_id_df_parasite$human_nonmito_homologs)] <- ""
gold_gene_accession_OG_id_df_parasite$present_human_nonmito <- as.numeric(gold_gene_accession_OG_id_df_parasite$human_nonmito_homologs != "")

# Get PFAM annotations
gold_gene_accession_OG_id_df_parasite_pfam <- merge(gold_gene_accession_OG_id_df_parasite, pfam_seprows, by.x = "gene_accession", by.y = "Accession", all.x = TRUE)
gold_gene_accession_OG_id_df_parasite_pfam <- gold_gene_accession_OG_id_df_parasite_pfam %>% group_by(gene_accession, taxid, PROG_id, present_Aca_mito, present_kinetoplastid_mito, present_Bdi_mito, present_n_parasite_clades, present_human_nonmito) %>% summarize(PfamList = paste0(sort(PfamList), collapse=","))
gold_gene_accession_OG_id_df_parasite_pfam$PfamList[is.na(gold_gene_accession_OG_id_df_parasite_pfam$PfamList)] <- ""

gold_gene_accession_OG_id_df_parasite_pfam$absense <- "Not powered"
gold_gene_accession_OG_id_df_parasite_pfam$absense[gold_gene_accession_OG_id_df_parasite_pfam$PROG_id %in% mito_PhROG_ids_detectable_human] <- "Powered"

# Get gene symbols
geneid_to_symbol <- read.delim(here("data", "orthogroups", "idmapping", "mitoepi_geneid_to_symbol_no_version.tsv"), header=FALSE)
colnames(geneid_to_symbol) <- c("protein_id", "symbol", "taxid")
gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_ids_notaxid <- gsub("^[0-9]+_", "", gold_gene_accession_OG_id_df_parasite_pfam$gene_accession)
gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_ids_notaxid <- gsub("\\.[0-9]{1}$", "", gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_ids_notaxid)
gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_taxid <- gsub("_.*", "", gold_gene_accession_OG_id_df_parasite_pfam$gene_accession)
gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_symbol <- geneid_to_symbol$symbol[match(gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_ids_notaxid, geneid_to_symbol$protein_id)]
gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_symbol[is.na(gold_gene_accession_OG_id_df_parasite_pfam$reference_protein_symbol)] <- ""

gold_gene_accession_OG_id_df_parasite_pfam <- gold_gene_accession_OG_id_df_parasite_pfam[order(gold_gene_accession_OG_id_df_parasite_pfam$present_n_parasite_clades, -gold_gene_accession_OG_id_df_parasite_pfam$present_human_nonmito, gold_gene_accession_OG_id_df_parasite_pfam$PROG_id, decreasing=TRUE),]

## Plot heatmap
gold_gene_accession_OG_id_df_parasite_pfam_select <- gold_gene_accession_OG_id_df_parasite_pfam %>% filter(taxid %in% selected_taxids)
gold_gene_accession_OG_id_df_parasite_pfam_select <- gold_gene_accession_OG_id_df_parasite_pfam_select %>% filter(present_n_parasite_clades >= 2) %>% filter(present_human_nonmito == 0) %>% filter(absense == "Powered")

# Build OGs matrix
species_tree_subtree <- keep.tip(species_tree, selected_taxids)
PROG_ids <- unique(gold_gene_accession_OG_id_df_parasite_pfam_select$PROG_id)
phrogs_long_eukaryota_parent_subset <- phrogs_long_eukaryota_parent %>% filter(taxid %in% species_tree_subtree$tip.label)
og_presence_mat <- matrix(data = 0, nrow = length(PROG_ids), ncol = length(species_tree_subtree$tip.label))
colnames(og_presence_mat) <- species_tree_subtree$tip.label
symbols <- c()
for (i in 1:length(PROG_ids)) {
  curr_OG_id <- PROG_ids[i]
  ogs_subset_long_curr <- phrogs_long_eukaryota_parent_subset %>% filter(PROG_id == curr_OG_id)
  ogs_subset_long_curr_nonmito <- ogs_subset_long_curr %>% filter(!protein_id %in% gold_gene_accession_OG_id_df$gene_accession)
  ogs_subset_long_curr_mito <- ogs_subset_long_curr %>% filter(protein_id %in% gold_gene_accession_OG_id_df$gene_accession)
  og_presence_mat[i, unique(ogs_subset_long_curr_nonmito$taxid)] <- 0.5
  og_presence_mat[i, unique(ogs_subset_long_curr_mito$taxid)] <- 1
  
  unique_pfams <- unique(unlist(strsplit(gold_gene_accession_OG_id_df_parasite_pfam_select$PfamList[gold_gene_accession_OG_id_df_parasite_pfam_select$PROG_id %in% curr_OG_id], ",")))
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
n_parasite_kingdoms <- as.numeric(rowSums(og_presence_mat_plot[,c("5689", "185431")]) == 2) + as.numeric(og_presence_mat_plot[,"1257118"] == 1) + as.numeric(og_presence_mat_plot[,"32595"] == 1)
og_presence_mat_reorder <- og_presence_mat_plot[order(og_presence_mat_plot[,"9606"], -n_parasite_kingdoms, -rowSums(og_presence_mat_plot)),]

my_palette <- colorRampPalette(c("white", "#999999", "#eb35ff"))(n = 100)
pdf("mito_PhROGs_parasite_targets_no.human.homolog_detectable_presence_mat.pdf", height=10, width=8)
heatmap.2(og_presence_mat_reorder, Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none", margins=c(8,16), na.rm=TRUE, col=my_palette, trace="none", key.xlab="Absence/Presence", cexRow = 0.4, cexCol = 1, key=FALSE, offsetRow=0)
dev.off()

# Venn diagram on a per PhROG basis
parasite_mito_set <- list('Kinetoplastida' = kinetoplast_OG_ids, 'Acanthamoeba' = mito_sets$Acanthamoeba, 'Babesia' = mito_sets$Babesia)
library(ggvenn)
pdf('parasite_PhROGs_venn.pdf', width=6, height=4)
ggvenn(parasite_mito_set, c("Kinetoplastida", "Acanthamoeba", "Babesia"), show_percentage = FALSE, set_name_size=4, text_size=6, fill_color = c("gray", "#00A36C", "purple"))
dev.off()

parasite_nohuman_powered_mito_set <- list('Kinetoplastida' = kinetoplast_OG_ids[!kinetoplast_OG_ids %in% phrogs_long_eukaryota_parent_human_collapse$PROG_id & kinetoplast_OG_ids %in% mito_PhROG_ids_detectable_human], 'Acanthamoeba' = mito_sets$Acanthamoeba[!mito_sets$Acanthamoeba %in% phrogs_long_eukaryota_parent_human_collapse$PROG_id & mito_sets$Acanthamoeba %in% mito_PhROG_ids_detectable_human], 'Babesia' = mito_sets$Babesia[!mito_sets$Babesia %in% phrogs_long_eukaryota_parent_human_collapse$PROG_id & mito_sets$Babesia %in% mito_PhROG_ids_detectable_human])
pdf('parasite_PhROGs_nohuman_powered_venn.pdf', width=6, height=4)
ggvenn(parasite_nohuman_powered_mito_set, c("Kinetoplastida", "Acanthamoeba", "Babesia"), show_percentage = FALSE, set_name_size=4, text_size=6, fill_color = c("gray", "#00A36C", "purple"))
dev.off()

## Write out
# write.table(gold_gene_accession_OG_id_df_parasite_pfam_select, here("data", "comparative", "mito_PhROGs_parasite_targets_pfam_absense.tsv"), sep="\t", row.names=FALSE)


