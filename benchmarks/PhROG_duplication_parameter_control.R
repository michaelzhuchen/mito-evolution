suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

# Set ggplot theme
theme_set(theme_classic())

# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
suffix <- args[1]

print(paste0("Analyzing ", suffix))

# Read in Eukaryota parent and child PhROGs
phrogs_long_eukaryota_parent <- read.table(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", suffix, "Node34_Eukaryota_parent_PhROGs_long.tsv"), sep="\t", header=FALSE)
colnames(phrogs_long_eukaryota_parent) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_parent_primary <- phrogs_long_eukaryota_parent %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id)) %>% ungroup()

phrogs_long_eukaryota <- read.table(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", suffix, "Node34_Eukaryota_PhROGs_long.tsv"), sep="\t", header=FALSE)
colnames(phrogs_long_eukaryota) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_primary <- phrogs_long_eukaryota %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id)) %>% ungroup()

# Filter for minimum of 4 species (same as LECA reconstruction)
n_minimum_species_leca <- 4
phrogs_long_eukaryota_parent_primary <- phrogs_long_eukaryota_parent_primary %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)
phrogs_long_eukaryota_primary <- phrogs_long_eukaryota_primary %>% group_by(PROG_id) %>% filter(length(unique(taxid)) >= n_minimum_species_leca)

# Get LECA PhROGs
phrogs_long_eukaryota_parent_primary_leca <- phrogs_long_eukaryota_parent_primary %>% filter(label == "Node34_Eukaryota")
eukaryota_parent_primary_leca_OG_protein_ids <- paste0(phrogs_long_eukaryota_parent_primary_leca$OG_id, "_", phrogs_long_eukaryota_parent_primary_leca$protein_id)
phrogs_long_eukaryota_primary_leca <- phrogs_long_eukaryota_primary %>% filter(paste0(OG_id, "_", protein_id) %in% eukaryota_parent_primary_leca_OG_protein_ids)

human_id2training <- read.table(here("data", "orthogroups", "idmapping", "human.id2training_uniprot_symbols.txt"), sep="\t", header=TRUE)
human_id2training <- human_id2training %>% mutate(protein_id = paste0("9606_", Entry))

# Read in Kay data and map using human gene ids
kay_data <- read.csv(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", "Kay_et_al_2026_Supplementary_Table_7.csv"), header=TRUE)
kay_data$LECA_cluster_id <- 1:nrow(kay_data)

# Collapse pathway annotation
selected_pathways <- c("replication.control", "Transcription", "DNA.repair", "Meiosis", "Nucleus.enclosure", "Nucleolus", "Lsm.Sm.complex.duplication", "Ubiquitin", "Proteasome", "endomembrane.development", "Translation", "Trafficking", "Splicing", "Mitochondrial.import")
indexes <- which(kay_data[,selected_pathways] == "Y", arr.ind = TRUE)
kay_data$pathway <- "Other"
kay_data$pathway[indexes[,1]] <- selected_pathways[indexes[,2]]

kay_data_human <- kay_data %>% group_by(HMM.or.Reference.protein) %>% filter() %>% ungroup() %>% filter(Node.type %in% "LECA") %>% filter(NCBI.gene.id != "-") %>% separate_rows(NCBI.gene.id, sep = " ")
kay_data_human <- kay_data_human %>% group_by(HMM.or.Reference.protein) %>% filter(all(NCBI.gene.id %in% human_id2training$HumanGeneID))
kay_data_human$protein_id <- human_id2training$protein_id[match(kay_data_human$NCBI.gene.id, human_id2training$HumanGeneID)]

kay_data_human$OG_PROG_id_parent_leca <- phrogs_long_eukaryota_parent_primary_leca$PROG_id[match(kay_data_human$protein_id, phrogs_long_eukaryota_parent_primary_leca$protein_id)]
kay_data_human$OG_PROG_id <- phrogs_long_eukaryota_primary$PROG_id[match(kay_data_human$protein_id, phrogs_long_eukaryota_primary$protein_id)]

# Bar plots by pathway
kay_data_human_leca_counts <- kay_data_human %>% group_by(HMM.or.Reference.protein) %>% filter(length(unique(LECA_cluster_id)) > 1) %>% mutate(n_LECA_genes = length(unique(OG_PROG_id[!is.na(OG_PROG_id)]))) %>% mutate(copies_label = case_when(n_LECA_genes > 1 ~ "Multicopy in LECA", n_LECA_genes == 1 ~ "Single copy in LECA", TRUE ~ "Not in LECA"))
kay_data_human_leca_counts <- kay_data_human_leca_counts %>% filter(copies_label %in% c("Single copy in LECA", "Multicopy in LECA"))
kay_data_human_leca_counts$copies_label <- factor(kay_data_human_leca_counts$copies_label, levels = c("Single copy in LECA", "Multicopy in LECA"))
kay_data_human_leca_counts_summary <- kay_data_human_leca_counts %>% group_by(pathway, copies_label) %>% summarize(count = n())
kay_data_human_leca_counts_pathway_order <- kay_data_human_leca_counts %>% group_by(pathway) %>% summarize(count = n()) %>% arrange(desc(count))
kay_data_human_leca_counts_summary$pathway <- factor(kay_data_human_leca_counts_summary$pathway, levels = kay_data_human_leca_counts_pathway_order$pathway)
pdf(paste0('barplot_Kay_2026_leca_copies_pathways_', suffix, '.pdf'), width = 6, height = 4)
ggplot(data = kay_data_human_leca_counts_summary, aes(x = pathway, y = count, fill = copies_label)) + geom_col() + coord_flip() + scale_fill_manual(values = c("gray", "black"))
dev.off()

# Bubble plot
leca_gene_counts <- kay_data_human %>% group_by(HMM.or.Reference.protein) %>% filter(length(unique(LECA_cluster_id)) > 1) %>% summarize(n_LECA_genes_kay = length(unique(LECA_cluster_id)), n_LECA_genes_chen = length(unique(OG_PROG_id[!is.na(OG_PROG_id)])))
leca_gene_counts_summary <- leca_gene_counts %>% group_by(n_LECA_genes_kay, n_LECA_genes_chen) %>% summarize(n_protein_families = length(unique(HMM.or.Reference.protein))) %>% ungroup()
pdf(paste0('bubbleplot_Kay_2026_leca_copies_', suffix, '.pdf'), width = 6, height = 4)
ggplot() + geom_point(data = leca_gene_counts_summary, aes(x = n_LECA_genes_kay, y = n_LECA_genes_chen, size = n_protein_families), shape=21, color="black", fill = "lightgray", alpha = 1) + scale_size_area(max_size = 10) + geom_line(data = data.frame(x = c(0,15), y = c(0,15)), aes(x = x, y = y), linetype = "dotted") + xlim(0,15) + ylim(0,15)
dev.off()

### OXPHOS negative control
## Map to human MitoPathways
mitopathways <- read.table(here("data", "annotation", "human_mitocarta", "human.path2gene.txt"), sep="\t", header=TRUE)
mitopathways_human_id2training <- merge(mitopathways, human_id2training, by.x="Gene", by.y="Symbol")
mitopathways_human_id2training$gene_accession <- paste0("9606_", mitopathways_human_id2training$Entry)
selected_mitopathways <- c("CI_subunits", "CII_subunits", "CIII_subunits", "CIV_subunits", "CV_subunits", "Mitochondrial_ribosome")
mitopathways_human_id2training_oxphos <- mitopathways_human_id2training %>% filter(MitoPathway %in% selected_mitopathways)
oxphos_protein_ids <- unique(mitopathways_human_id2training_oxphos$protein_id)

# Filter by OXPHOS Eukaryota parent PhROGs to specifically find duplications
leca_oxphos_ogs <- phrogs_long_eukaryota_primary %>% group_by(OG_id) %>% filter(any(protein_id %in% oxphos_protein_ids)) %>% ungroup()
phrogs_long_eukaryota_parent_primary_leca_oxphos <- phrogs_long_eukaryota_parent_primary_leca %>% group_by(OG_id, PROG_id) %>% filter(any(protein_id %in% oxphos_protein_ids)) %>% ungroup()
leca_oxphos_ogs$PROG_id[!leca_oxphos_ogs$protein_id %in% phrogs_long_eukaryota_parent_primary_leca_oxphos$protein_id] <- NA

leca_oxphos_ogs$human_symbol <- human_id2training$Symbol[match(leca_oxphos_ogs$protein_id, human_id2training$protein_id)]
leca_oxphos_ogs$pathway <- mitopathways_human_id2training_oxphos$MitoPathway[match(leca_oxphos_ogs$protein_id, mitopathways_human_id2training_oxphos$protein_id)]
leca_oxphos_ogs_counts <- leca_oxphos_ogs %>% group_by(OG_id) %>% mutate(n_LECA_genes = length(unique(PROG_id[!is.na(PROG_id)])), n_species_any_copy = length(unique(taxid)), n_species_multiple_copies = length(unique(taxid[duplicated(taxid)]))) %>% mutate(fraction_species_multiple_copies = n_species_multiple_copies / n_species_any_copy) %>% group_by(OG_id, n_LECA_genes, n_species_any_copy, n_species_multiple_copies, fraction_species_multiple_copies) %>% filter(protein_id %in% oxphos_protein_ids) %>% summarize(human_protein_ids = paste0(protein_id, collapse=","), human_symbols = paste0(sort(human_symbol), collapse=","), pathway = paste0(unique(pathway), collapse=",")) %>% separate_rows(human_symbols, sep=",")

leca_oxphos_ogs_counts$label <- "Not in LECA"
leca_oxphos_ogs_counts$label[leca_oxphos_ogs_counts$n_LECA_genes == 1] <- "Single copy in LECA"
leca_oxphos_ogs_counts$label[leca_oxphos_ogs_counts$n_LECA_genes > 1] <- "Multi copy in LECA"
leca_oxphos_ogs_counts <- leca_oxphos_ogs_counts %>% filter(label %in% c("Single copy in LECA", "Multi copy in LECA"))
leca_oxphos_ogs_counts$pathway <- factor(leca_oxphos_ogs_counts$pathway, levels = rev(selected_mitopathways))
pdf(paste0("barplot_oxphos_leca_copies_", suffix, ".pdf"), width = 6, height = 4)
ggplot(data = leca_oxphos_ogs_counts, aes(x = pathway, fill = label)) + geom_bar() + coord_flip() + scale_fill_manual(values = c("black", "gray"))
dev.off()


# Print stats
print((paste0("Stats for ", suffix)))
sum(kay_data_human_leca_counts$copies_label == "Multicopy in LECA") / nrow(kay_data_human_leca_counts)
length(unique(kay_data_human_leca_counts$HMM.or.Reference.protein[kay_data_human_leca_counts$copies_label == "Multicopy in LECA"])) / length(unique(kay_data_human_leca_counts$HMM.or.Reference.protein))
length(unique(leca_gene_counts$HMM.or.Reference.protein[leca_gene_counts$n_LECA_genes_chen == leca_gene_counts$n_LECA_genes_kay])) / length(unique(leca_gene_counts$HMM.or.Reference.protein))
cor(leca_gene_counts$n_LECA_genes_kay, leca_gene_counts$n_LECA_genes_chen, method = "spearman")
sum(leca_oxphos_ogs_counts$label == "Single copy in LECA") / nrow(leca_oxphos_ogs_counts)
length(unique(leca_oxphos_ogs_counts$OG_id[leca_oxphos_ogs_counts$label == "Single copy in LECA"])) / length(unique(leca_oxphos_ogs_counts$OG_id))


## Plot duplication stats across species overlap thresholds
duplication_stats <- read.csv(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", "duplication_control_stats.csv"))
pdf("scatterplot_duplication_stats.pdf", width = 4, height = 4)
ggplot(data = duplication_stats, aes(x = fraction.Kay.families.with.exactly.correct.copies.in.LECA, y = fraction.OXPHOS.single.copy.families.in.LECA, label = species.overlap.threshold)) + geom_point(size=4) + geom_text_repel(size=4) + xlim(0.6,1) + ylim(0.6,1) + xlab("Recall (fraction of duplicated LECA families with same copies as literature)") + ylab("Precision (fraction of OxPhos LECA families with single copy)") + coord_fixed(ratio = 1)
dev.off()
