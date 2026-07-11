library(tidyverse)

virus_taxids <- c("212035", "554168", "2080449")

# Read in mapping for ACANB and ACABI ids
ACANB_best_transcript_isoforms <- read.delim(here("data", "orthogroups", "idmapping", "ACANBgene.genes_best.isoform.txt"))$transcriptID

# Read in mimivirus hmmer results
tophit_mimivirus <- read.table(here("specific_analyses_per_organism", "acanthamoeba_viral_homologs", "OG_hmmsearch_mimivirus_expect1e-3_DBSIZE.379668_tophit.tsv"))
colnames(tophit_mimivirus) <- c("accession", "Orthogroup")
tophit_mimivirus <- tophit_mimivirus %>% mutate(taxid = gsub("_.*", "", accession))
tophit_mimivirus$BOOL_PRIMARY_OG <- TRUE

combined_fusion_protein_filter_disjoint_mimivirus <- read.table(here("specific_analyses_per_organism", "acanthamoeba_viral_homologs", "OG_hmmsearch_mimivirus_expect1e-3_DBSIZE.379668_fusions.tsv"), sep="\t", header=FALSE)
colnames(combined_fusion_protein_filter_disjoint_mimivirus) <- c("accession", "Orthogroup")
combined_fusion_protein_filter_disjoint_mimivirus <- combined_fusion_protein_filter_disjoint_mimivirus %>% separate_rows(Orthogroup, sep=", ")
combined_fusion_protein_filter_disjoint_mimivirus <- combined_fusion_protein_filter_disjoint_mimivirus %>% mutate(taxid = gsub("_.*", "", accession))
combined_fusion_protein_filter_disjoint_mimivirus$BOOL_PRIMARY_OG <- FALSE

mimivirus_ogs <- rbind(tophit_mimivirus, combined_fusion_protein_filter_disjoint_mimivirus)
mimivirus_ogs <- mimivirus_ogs %>% distinct(accession, Orthogroup, .keep_all = TRUE)

# Read in OGs
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")
ogs_long_intersect_aca_mimivirus <- ogs_long %>% filter(Orthogroup %in% mimivirus_ogs$Orthogroup) %>% filter(taxid == "1257118")
intersect_aca_mimivirus_og_ids <- unique(ogs_long_intersect_aca_mimivirus$Orthogroup)

ogs_long_select <- ogs_long %>% filter(Orthogroup %in% intersect_aca_mimivirus_og_ids)
mimivirus_ogs_select <- mimivirus_ogs %>% filter(Orthogroup %in% intersect_aca_mimivirus_og_ids)
ogs_long_mimivirus <- rbind(ogs_long_select, mimivirus_ogs_select)

# Pie chart of ACANB genes with/without viral homolog
ogs_long_select_acanb <- ogs_long_select %>% filter(taxid %in% c("1257118"))
n_acanb_genes <- length(ACANB_best_transcript_isoforms)
n_acanb_genes_mimivirus_homologs <- length(unique(ogs_long_select_acanb$accession))
n_acanb_genes_no_mimivirus_homologs <- n_acanb_genes - n_acanb_genes_mimivirus_homologs
pdf('pie_mimivirus_ACANB_homologs.pdf', width = 10, height = 4)
pie(c(n_acanb_genes_no_mimivirus_homologs, n_acanb_genes_mimivirus_homologs), labels = c(paste0("No viral orthogroup homolog", " (", n_acanb_genes_no_mimivirus_homologs, ")"), paste0("Viral orthogroup homolog", " (", n_acanb_genes_mimivirus_homologs, ")")))
dev.off()

# Pie chart of ACANB genes viral sister
mimivirus_to_aca <- read.delim(here("specific_analyses_per_organism", "acanthamoeba_viral_homologs", "HGT_results_virus_to_1257118.tsv"), sep="\t", header=TRUE)
mimivirus_to_aca <- mimivirus_to_aca %>% mutate(ID = gsub("1257118_", "", recipient_protein))
mimivirus_to_aca_sister <- mimivirus_to_aca %>% filter(selected_donor_sister_tiplabels != "")
ogs_long_select_acanb <- ogs_long_select %>% filter(taxid %in% c("1257118"))
n_acanb_genes_mimivirus_sister <- length(unique(mimivirus_to_aca_sister$ID))
n_acanb_genes_no_mimivirus_sister <- n_acanb_genes_mimivirus_homologs - n_acanb_genes_mimivirus_sister
pdf('pie_mimivirus_ACANB_sister.pdf', width = 6, height = 4)
pie(c(n_acanb_genes_no_mimivirus_sister, n_acanb_genes_mimivirus_sister), labels = c(paste0("No viral sister homolog", " (", n_acanb_genes_no_mimivirus_sister, ")"), paste0("Viral sister homolog", " (", n_acanb_genes_mimivirus_sister, ")")))
dev.off()

## Write out results table
acanb_homologs <- data.frame(ID = ACANB_best_transcript_isoforms)
acanb_homologs <- acanb_homologs %>% mutate(Protein_ID = paste0("1257118_", ID))

ogs_long_mimivirus_collapse <- ogs_long_mimivirus %>% group_by(Orthogroup) %>% mutate(viral_orthogroup_homologs = paste0(sort(unique(accession[taxid %in% virus_taxids])), collapse = ",")) %>% filter(taxid == "1257118") %>% group_by(accession) %>% summarize(viral_orthogroup_homologs = paste0(sort(unique(unlist(strsplit(viral_orthogroup_homologs, split=",")))), collapse=","))
acanb_homologs_mimivirus <- merge(acanb_homologs, ogs_long_mimivirus_collapse, by.x="Protein_ID", by.y="accession", all.x=TRUE, all.y=FALSE)
acanb_homologs_mimivirus_sister <- merge(acanb_homologs_mimivirus, mimivirus_to_aca_sister, by.x="Protein_ID", by.y="recipient_protein", all.x=TRUE, all.y=FALSE)

acanb_homologs_mimivirus_sister_reduce <- acanb_homologs_mimivirus_sister[,c("Protein_ID", "viral_orthogroup_homologs", "selected_donor_sister_tiplabels", "sister_node_support")]
colnames(acanb_homologs_mimivirus_sister_reduce) <- c("Protein_ID", "Viral_homologs", "Viral_homologs_in_sister_clade", "Bootstrap_support")
write.csv(acanb_homologs_mimivirus_sister_reduce, here("ACANB_mimivirus_homologs.csv"), row.names = FALSE)

