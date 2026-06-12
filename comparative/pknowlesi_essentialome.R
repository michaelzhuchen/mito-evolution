library(tidyverse)

# Read in P. knowlesi data
pkno_essential <- read.csv(here("data", "annotation", "essentialome", "PknoH_TPNessentiality&OG(orthology)_030725.csv"))
pkno_id_mapping <- read.delim(here("data", "orthogroups", "idmapping", "pkno_genename.to.uniprot_idmapping_2025_03_13.tsv"), sep="\t", header=TRUE)
pkno_essential_uniprot <- merge(pkno_essential, pkno_id_mapping, by.x="Gene.ID", by.y="From")
pkno_essential_uniprot$P.knowlesi.strain.H.piggyBac.mutagenesis...Hybrid.model.score[pkno_essential_uniprot$P.knowlesi.strain.H.piggyBac.mutagenesis...Hybrid.model.score == "N/A"] <- NA
pkno_essential_uniprot$P.knowlesi.strain.H.piggyBac.mutagenesis...Hybrid.model.score <- as.numeric(pkno_essential_uniprot$P.knowlesi.strain.H.piggyBac.mutagenesis...Hybrid.model.score)

gold_gene_accession_OG_id_df <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_2026.04.05.tsv"), sep="\t", header=TRUE)

phrogs_long_eukaryota_parent <- read.table(here("data", "phylogenetically_resolved_orthogroups", "species_tree_1", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long_eukaryota_parent) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
phrogs_long_eukaryota_parent <- phrogs_long_eukaryota_parent %>% filter(BOOL_primary_OG) %>% mutate(OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id)) %>% ungroup()

# Essentiality with PhROGs
phrogs_long_eukaryota_parent_pkno <- phrogs_long_eukaryota_parent %>% filter(taxid == "5851")

# Read in OGs
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")

# Add OGs for proteins not in PhROGs
ogs_long_pkno <- ogs_long %>% filter(taxid == "5851") %>% filter(BOOL_PRIMARY_OG) %>% filter(!accession %in% phrogs_long_eukaryota_parent_pkno$protein_id)
ogs_long_pkno$PROG_id <- ogs_long_pkno$Orthogroup
ogs_long_pkno$protein_id <- ogs_long_pkno$accession
phrogs_long_eukaryota_parent_pkno <- rbind(phrogs_long_eukaryota_parent_pkno[,c("PROG_id", "protein_id")], ogs_long_pkno[,c("PROG_id", "protein_id")])

phrogs_long_eukaryota_parent_pkno$pkno_hybrid_model_score <- pkno_essential_uniprot$P.knowlesi.strain.H.piggyBac.mutagenesis...Hybrid.model.score[match(gsub("5851_", "", phrogs_long_eukaryota_parent_pkno$protein_id), pkno_essential_uniprot$Entry.Name)]

phrogs_long_eukaryota_parent_pkno <- phrogs_long_eukaryota_parent_pkno %>% filter(!is.na(pkno_hybrid_model_score))

completed_mitoproteomes_species_list <- c("9606", "559292", "3702", "1257118", "5689", "185431", "5741", "32595")

phrogs_long_eukaryota_parent_all <- phrogs_long_eukaryota_parent %>% filter(taxid %in% completed_mitoproteomes_species_list)
all_PhROG_ids <- unique(phrogs_long_eukaryota_parent_all$PROG_id)

gold_gene_accession_OG_id_subset_df <- gold_gene_accession_OG_id_df[gold_gene_accession_OG_id_df$taxid %in% completed_mitoproteomes_species_list,]
gold_gene_accession_OG_id_subset_df %>% filter(OG_id != "") %>% group_by(taxid) %>% summarize(n=n())
gold_gene_accession_OG_id_subset_df <- gold_gene_accession_OG_id_subset_df %>% filter(OG_id != "")
phrogs_long_eukaryota_parent_mito <- phrogs_long_eukaryota_parent %>% filter(protein_id %in% gold_gene_accession_OG_id_subset_df$gene_accession)
mito_PhROG_ids <- unique(phrogs_long_eukaryota_parent_mito$PROG_id)

phrogs_long_eukaryota_parent_pkno$category <- "Nonmito"
phrogs_long_eukaryota_parent_pkno$category[phrogs_long_eukaryota_parent_pkno$PROG_id %in% mito_PhROG_ids] <- "Mito"

phrogs_long_eukaryota_parent_pkno_summary <- phrogs_long_eukaryota_parent_pkno %>% filter(PROG_id %in% all_PhROG_ids) %>% group_by(PROG_id) %>% summarize(pkno_hybrid_model_score = mean(pkno_hybrid_model_score, na.rm = TRUE), pkno_genes = paste0(protein_id[grepl("5851_", protein_id, fixed=TRUE)], collapse = ","), category = unique(category))

pkno_essential_threshold <- 0.26 # essential are less than this value
pkno_dispensible_threshold <- 0.88 # dispensible are greater than this value

# Add mito and parasite drug targets
gold_gene_accession_OG_id_df_parasite_pfam_select <- read.table(here("data", "comparative", "mito_PhROGs_parasite_targets_pfam_absense.tsv"), header=TRUE)
phrogs_long_eukaryota_parent_pkno_summary$label <- "Other non-mito"
phrogs_long_eukaryota_parent_pkno_summary$label[which(phrogs_long_eukaryota_parent_pkno_summary$PROG_id %in% mito_PhROG_ids)] <- "Other mito"
phrogs_long_eukaryota_parent_pkno_summary$label[which(phrogs_long_eukaryota_parent_pkno_summary$PROG_id %in% unique(gold_gene_accession_OG_id_df_parasite_pfam_select$PROG_id))] <- "Candidate target"
phrogs_long_eukaryota_parent_pkno_summary$label <- factor(phrogs_long_eukaryota_parent_pkno_summary$label, levels = c("Other non-mito", "Other mito", "Candidate target"))

pdf('pkno_essential_mito_drug_candidates_PhROGs.pdf', height=2, width=5)
ggplot(data = phrogs_long_eukaryota_parent_pkno_summary, aes(x = pkno_hybrid_model_score, color = label)) + geom_density(bw = 0.1) + xlab("Plasmodium knowlesi hybrid model score") + geom_vline(aes(xintercept = pkno_essential_threshold), linetype = "dashed") + geom_vline(aes(xintercept = pkno_dispensible_threshold), linetype = "dashed") + scale_color_manual(values=c("gray", "#0047AB", "#eb35ff"))
dev.off()
