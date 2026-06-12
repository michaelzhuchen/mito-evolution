
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(gplots)
library(basetheme)
library(RColorBrewer)
library(castor)
library(ape)
library(phytools)

# Set ggplot theme
theme_set(theme_classic())

# Read in species tree
species_tree <- read.tree(here("data/species_phylogeny/processed_species_tree", paste0("species_tree_1.nwk")))
species_tree_labels <- c(species_tree$tip.label, species_tree$node.label)
archaeplastida_subtree <- get_subtree_at_node(species_tree, "Node39_Archaeplastida")$subtree
archaeplastida_subtree_labels <- c(archaeplastida_subtree$tip.label, archaeplastida_subtree$node.label)

progs_long_agg <- c()

# Read in PhROGs
progs_long <- read.table(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", "HGT.nomaxspp", "Node34_Eukaryota_parent_PhROGs_long.tsv"), sep="\t", header=FALSE)
colnames(progs_long) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
progs_long$method <- "No max recipients"
progs_long_agg <- rbind(progs_long_agg, progs_long)
progs_long <- read.table(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", "HGT.100maxspp", "Node34_Eukaryota_parent_PhROGs_long.tsv"), sep="\t", header=FALSE)
colnames(progs_long) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
progs_long$method <- "Max 100 recipients"
progs_long_agg <- rbind(progs_long_agg, progs_long)
progs_long <- read.table(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", "HGT.40maxspp", "Node34_Eukaryota_parent_PhROGs_long.tsv"), sep="\t", header=FALSE)
colnames(progs_long) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
progs_long$method <- "Max 40 recipients"
progs_long_agg <- rbind(progs_long_agg, progs_long)
progs_long <- read.table(here("data", "benchmarks", "phylogenetically_resolved_orthogroups", "HGT.10maxspp", "Node34_Eukaryota_parent_PhROGs_long.tsv"), sep="\t", header=FALSE)
colnames(progs_long) <- c("PROG_id", "protein_id", "label", "mito_localization_prob", "BOOL_NONVERTICAL", "BOOL_primary_OG")
progs_long$method <- "Max 10 recipients"
progs_long_agg <- rbind(progs_long_agg, progs_long)

progs_long <- progs_long_agg %>% mutate(OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% filter(BOOL_primary_OG) %>% group_by(method, OG_id) %>% filter(!duplicated(protein_id))

## Read in files for uniprot id mapping
ath_idmapping <- read.delim(here("data", "orthogroups", "idmapping", "ath_idmapping_2024_11_08.tsv"))
orthogroups_collapse_accessions_seprows_uniprot <- read.table(here("data", "orthogroups", "idmapping", "refined_OGs_euk673spp_uniprot_mapping_wide.txt"), sep="\t", header=FALSE)
colnames(orthogroups_collapse_accessions_seprows_uniprot) <- c("OG_id", "accessions", "uniprot_accession")
# Read in new Tbr mapping for fasta with best isoform per gene
tbr_mapping <- read.table(here("data", "orthogroups", "idmapping", "tbrgene.mapping.txt"), header=TRUE)
tbr_mapping$transcriptID <- gsub(":", "_", tbr_mapping$transcriptID, fixed=TRUE)
tbr_mapping$ID <- paste0("185431_", tbr_mapping$ID)
tbr_mapping$transcriptID <- paste0("185431_", tbr_mapping$transcriptID)
orthogroups_collapse_accessions_seprows_uniprot$accessions[orthogroups_collapse_accessions_seprows_uniprot$accessions %in% tbr_mapping$transcriptID] <- tbr_mapping$ID[match(orthogroups_collapse_accessions_seprows_uniprot$accessions[orthogroups_collapse_accessions_seprows_uniprot$accessions %in% tbr_mapping$transcriptID], tbr_mapping$transcriptID)]


mitoribo_table <- read.csv(here("data", "annotation", "pathways", "mitoribosome_structure_annotation.csv"), header=TRUE)
mitoribo_table_core_bacterial <- mitoribo_table[which(mitoribo_table$type == "core bacterial" & mitoribo_table$sum == 4),]
mitoribo_table_shared_euk <- mitoribo_table[which(mitoribo_table$type == "shared euk" & mitoribo_table$sum == 4),]
mitoribo_table_core <- rbind(mitoribo_table_core_bacterial, mitoribo_table_shared_euk)
ath_mitoribo <- mitoribo_table_core$Arabidopsis
ath_mitoribo <- ath_mitoribo[ath_mitoribo != ""]
# Retrieve all possible Arabidopsis IDs, since some annotated IDs are not the canonical entry
ath_geneid <- ath_idmapping$From[match(ath_mitoribo, ath_idmapping$Entry)]
ath_uniprot_ids <- ath_idmapping$Entry[match(ath_geneid, ath_idmapping$From)]
ath_mitoribo_accessions <- orthogroups_collapse_accessions_seprows_uniprot$accessions[match(ath_uniprot_ids, orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession)]

spinach_chlororibo <- mitoribo_table$Spinacea_chloroplast_ribosome_6ERI
spinach_chlororibo <- spinach_chlororibo[spinach_chlororibo != ""]
chlororibo_accessions <- orthogroups_collapse_accessions_seprows_uniprot$accessions[match(spinach_chlororibo, orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession)]

ps1_table <- read.csv(here("data", "annotation", "pathways", "PSI_structure_annotation.csv"), header=TRUE)
ps1_table <- ps1_table %>% filter(type == "core")
ps1_accessions <- orthogroups_collapse_accessions_seprows_uniprot$accessions[match(ps1_table$uniprot_accession, orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession)]
ps2_table <- read.csv(here("data", "annotation", "pathways", "PSII_structure_annotation.csv"), header=TRUE)
ps2_table <- ps2_table %>% filter(type == "core")
ps2_accessions <- orthogroups_collapse_accessions_seprows_uniprot$accessions[match(ps2_table$uniprot_accession, orthogroups_collapse_accessions_seprows_uniprot$uniprot_accession)]

# Select protein sets
plastid_proteins <- unique(c(ath_mitoribo_accessions, chlororibo_accessions, ps1_accessions, ps2_accessions))

progs_long_plastid <- progs_long %>% filter(protein_id %in% plastid_proteins)
progs_long_plastid$complex <- ""
progs_long_plastid$complex[progs_long_plastid$protein_id %in% ath_mitoribo_accessions] <- "Core mitochondrial ribosome (48)"
progs_long_plastid$complex[progs_long_plastid$protein_id %in% chlororibo_accessions] <- "Chloroplast ribosome (53)"
progs_long_plastid$complex[progs_long_plastid$protein_id %in% ps1_accessions] <- "Photosystem I (13)"
progs_long_plastid$complex[progs_long_plastid$protein_id %in% ps2_accessions] <- "Photosystem II (16)"
progs_long_plastid_summary <- progs_long_plastid %>% group_by(method, complex, label) %>% summarize(n = n(), .groups = "drop")

# Collapse labels and get counts
progs_long_plastid_summary$label_group <- progs_long_plastid_summary$label
progs_long_plastid_summary$label_group[progs_long_plastid_summary$label %in% c("Node39_Archaeplastida", "Node40_Viridiplantae", "Node41_Streptophyta", "Node42_Streptophytina", "Node43_Streptophytina")] <- "Within ancestral plants"
progs_long_plastid_summary$label_group[!progs_long_plastid_summary$label %in% c("Node34_Eukaryota", "Node35_Diphoda", "Node36_Diaphorectickes", "Node37_CAM_Haptista", "Node38_CAM", "Node39_Archaeplastida", "Node40_Viridiplantae", "Node41_Streptophyta", "Node42_Streptophytina", "Node43_Streptophytina", "Within ancestral plants")] <- "Within land plants"
progs_long_plastid_summary_group <- progs_long_plastid_summary %>% group_by(method, complex, label_group) %>% summarize(n_proteins = sum(n), .groups = "drop")

# Naive LCA on PhROGs
progs_long_plastid_LCA <- progs_long %>% filter(method == "Max 10 recipients") %>% filter(PROG_id %in% progs_long_plastid$PROG_id) %>% group_by(PROG_id) %>% mutate(mrca_index = get_mrca_of_set(species_tree, unique(taxid))) %>% mutate(label = species_tree_labels[mrca_index]) %>% filter(protein_id %in% plastid_proteins)
progs_long_plastid_LCA$complex <- ""
progs_long_plastid_LCA$complex[progs_long_plastid_LCA$protein_id %in% ath_mitoribo_accessions] <- "Core mitochondrial ribosome (48)"
progs_long_plastid_LCA$complex[progs_long_plastid_LCA$protein_id %in% chlororibo_accessions] <- "Chloroplast ribosome (53)"
progs_long_plastid_LCA$complex[progs_long_plastid_LCA$protein_id %in% ps1_accessions] <- "Photosystem I (13)"
progs_long_plastid_LCA$complex[progs_long_plastid_LCA$protein_id %in% ps2_accessions] <- "Photosystem II (16)"
progs_long_plastid_LCA$label_group <- progs_long_plastid_LCA$label
progs_long_plastid_LCA$label_group[progs_long_plastid_LCA$label %in% c("Node39_Archaeplastida", "Node40_Viridiplantae", "Node41_Streptophyta", "Node42_Streptophytina", "Node43_Streptophytina")] <- "Within ancestral plants"
progs_long_plastid_LCA$label_group[!progs_long_plastid_LCA$label %in% c("Node34_Eukaryota", "Node35_Diphoda", "Node36_Diaphorectickes", "Node37_CAM_Haptista", "Node38_CAM", "Node39_Archaeplastida", "Node40_Viridiplantae", "Node41_Streptophyta", "Node42_Streptophytina", "Node43_Streptophytina", "Within ancestral plants")] <- "Within land plants"
progs_long_plastid_LCA_summary_group <- progs_long_plastid_LCA %>% group_by(complex, label_group) %>% summarize(n_proteins = n(), .groups = "drop")

# Naive LCA on OGs
ogs_long_plastid_LCA <- progs_long %>% filter(method == "Max 10 recipients") %>% filter(OG_id %in% progs_long_plastid$OG_id) %>% group_by(OG_id) %>% mutate(mrca_index = get_mrca_of_set(species_tree, unique(taxid))) %>% mutate(label = species_tree_labels[mrca_index]) %>% filter(protein_id %in% plastid_proteins)
ogs_long_plastid_LCA$complex <- ""
ogs_long_plastid_LCA$complex[ogs_long_plastid_LCA$protein_id %in% ath_mitoribo_accessions] <- "Core mitochondrial ribosome (48)"
ogs_long_plastid_LCA$complex[ogs_long_plastid_LCA$protein_id %in% chlororibo_accessions] <- "Chloroplast ribosome (53)"
ogs_long_plastid_LCA$complex[ogs_long_plastid_LCA$protein_id %in% ps1_accessions] <- "Photosystem I (13)"
ogs_long_plastid_LCA$complex[ogs_long_plastid_LCA$protein_id %in% ps2_accessions] <- "Photosystem II (16)"
ogs_long_plastid_LCA$label_group <- ogs_long_plastid_LCA$label
ogs_long_plastid_LCA$label_group[ogs_long_plastid_LCA$label %in% c("Node39_Archaeplastida", "Node40_Viridiplantae", "Node41_Streptophyta", "Node42_Streptophytina", "Node43_Streptophytina")] <- "Within ancestral plants"
ogs_long_plastid_LCA$label_group[!ogs_long_plastid_LCA$label %in% c("Node34_Eukaryota", "Node35_Diphoda", "Node36_Diaphorectickes", "Node37_CAM_Haptista", "Node38_CAM", "Node39_Archaeplastida", "Node40_Viridiplantae", "Node41_Streptophyta", "Node42_Streptophytina", "Node43_Streptophytina", "Within ancestral plants")] <- "Within land plants"
ogs_long_plastid_LCA_summary_group <- ogs_long_plastid_LCA %>% group_by(complex, label_group) %>% summarize(n_proteins = n(), .groups = "drop")

# Combine
progs_long_plastid_LCA_summary_group$method <- "PhROG (prok->euk HGT-aware)"
ogs_long_plastid_LCA_summary_group$method <- "OG (HGT-naive)"
progs_long_plastid_summary_combined <- rbind(progs_long_plastid_summary_group, progs_long_plastid_LCA_summary_group, ogs_long_plastid_LCA_summary_group)

progs_long_plastid_summary_combined$complex <- factor(progs_long_plastid_summary_combined$complex, levels = c("Core mitochondrial ribosome (48)", "Chloroplast ribosome (53)", "Photosystem I (13)", "Photosystem II (16)"))
progs_long_plastid_summary_combined$method <- factor(progs_long_plastid_summary_combined$method, levels = c("No max recipients", "Max 100 recipients", "Max 40 recipients", "Max 10 recipients", "PhROG (prok->euk HGT-aware)", "OG (HGT-naive)"))

progs_long_plastid_summary_combined$label_group <- factor(progs_long_plastid_summary_combined$label_group, levels = sort(unique(progs_long_plastid_summary_combined$label_group), decreasing=TRUE))

progs_long_plastid_summary_combined_nomissing <- rbind(progs_long_plastid_summary_combined, cbind(expand.grid(method=levels(progs_long_plastid_summary_combined$method), complex=levels(progs_long_plastid_summary_combined$complex), label_group=levels(progs_long_plastid_summary_combined$label_group), n_proteins=0)))
progs_long_plastid_summary_combined_nomissing <- progs_long_plastid_summary_combined_nomissing %>% distinct(method, complex, label_group, .keep_all = TRUE)

# Plot
pdf('plastid_inferred_gain_timing_by_complex.pdf', width = 8, height = 6)
ggplot(data = progs_long_plastid_summary_combined_nomissing, aes(x = label_group, y = n_proteins, fill = method)) + geom_col(position = "dodge", width = 0.7) + scale_fill_manual(values = c("#00A36C40", "#00A36C60", "#00A36C80", "#00A36C", "gray", "lightgray")) + coord_flip() + facet_wrap(~ complex, nrow = 1) + xlab("Inferred acquisition timing") + ylab("Number of proteins") + theme(legend.position = "bottom")
dev.off()
