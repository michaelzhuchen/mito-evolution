
library(tidyverse)
theme_set(theme_classic())


# Read in OGs
ogs_long <- read.table(here("data", "orthogroups", "refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")

# Get PFAM annotations
pfam_seprows <- read.table(here("data", "annotation", "pfam", "euk673spp_prokgroups_pfam_list.filtered_seprows.tsv"), sep="\t", header=TRUE)

## Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data", "taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
euk_taxids <- uniprot_proteomes_tax$TaxId[which(uniprot_proteomes_tax$domain == "Eukaryota")]

# Select OGs to resample
OG_ids <- c("OG0005029", "OG0011422")

# Identify domains to resample
selected_pfam_domains_threshold <- 0.3

ogs_long_select <- ogs_long %>% filter(Orthogroup %in% OG_ids)

ogs_long_select_pfam <- merge(ogs_long_select, pfam_seprows, by.x="accession", by.y="Accession")

OG_domains <- c()

for (curr_OG_id in OG_ids) {
  
  n_proteins <- length(unique(ogs_long_select$accession[ogs_long_select$Orthogroup == curr_OG_id]))
  
  ogs_long_select_pfam_summary <- ogs_long_select_pfam %>% filter(Orthogroup == curr_OG_id) %>% group_by(PfamList) %>% summarize(fraction_proteins_with_domain = length(unique(accession)) / n_proteins, .groups = "drop")
  
  ogs_long_select_pfam_summary <- ogs_long_select_pfam_summary[order(ogs_long_select_pfam_summary$fraction_proteins_with_domain, decreasing=TRUE),]
  selected_pfam_domains <- ogs_long_select_pfam_summary$PfamList[ogs_long_select_pfam_summary$fraction_proteins_with_domain > selected_pfam_domains_threshold]
  
  ogs_long_select_pfam_architecture <- ogs_long_select_pfam %>% filter(Orthogroup == curr_OG_id) %>% filter(PfamList %in% selected_pfam_domains) %>% group_by(Orthogroup, accession) %>% filter(!duplicated(PfamList)) %>% summarize(PfamList_collapse = paste0(selected_pfam_domains[selected_pfam_domains %in% PfamList], collapse = ","), .groups = "drop")
  ogs_long_select_pfam_architecture_summary <- ogs_long_select_pfam_architecture %>% group_by(Orthogroup, PfamList_collapse) %>% summarize(fraction_with_domains = length(unique(accession)) / n_proteins, .groups = "drop") %>% rowwise() %>% mutate(weighted_fraction_with_domains = fraction_with_domains * length(unlist(strsplit(PfamList_collapse, split=","))))
  ogs_long_select_pfam_architecture_summary <- ogs_long_select_pfam_architecture_summary[order(ogs_long_select_pfam_architecture_summary$weighted_fraction_with_domains, decreasing=TRUE),]
  
  ogs_long_select_pfam_architecture_summary$rank <- 1:nrow(ogs_long_select_pfam_architecture_summary)
  ogs_long_select_pfam_architecture_summary_seprows <- ogs_long_select_pfam_architecture_summary %>% separate_rows(PfamList_collapse, sep=",")
  ogs_long_select_pfam_architecture_summary_seprows_nodups <- ogs_long_select_pfam_architecture_summary_seprows %>% filter(!duplicated(PfamList_collapse))
  
  selected_ranks <- unique(ogs_long_select_pfam_architecture_summary_seprows_nodups$rank)
  
  for (rank in selected_ranks) {
    OG_domain <- data.frame(OG_id = curr_OG_id, new_OG_id = paste0(curr_OG_id, "_", rank), selected_pfam_domains = ogs_long_select_pfam_architecture_summary$PfamList_collapse[ogs_long_select_pfam_architecture_summary$rank == rank])
    OG_domains <- rbind(OG_domains, OG_domain)
  }
}


ogs_long_to_split <- ogs_long %>% filter(Orthogroup %in% OG_domains$OG_id)
ogs_long_to_split_pfam <- merge(ogs_long_to_split, pfam_seprows, by.x="accession", by.y="Accession")

resampled_OGs_long <- c()

# Retrieve resampled OGs by consensus pfam
for (i in 1:nrow(OG_domains)) {
  OG_domains_curr <- OG_domains[i,]
  
  ogs_long_select <- ogs_long_to_split %>% filter(Orthogroup %in% OG_domains_curr$OG_id)
  ogs_long_select_pfam <- ogs_long_to_split_pfam %>% filter(Orthogroup %in% OG_domains_curr$OG_id)
  
  ogs_long_select_pfam_summary <- ogs_long_select_pfam %>% group_by(PfamList) %>% summarize(fraction_proteins_with_domain = length(unique(accession)) / length(unique(ogs_long_select$accession)))
  
  selected_pfam_domains <- unlist(strsplit(OG_domains_curr$selected_pfam, split=","))
  
  ogs_long_select_pfam_has_domains <- ogs_long_select_pfam %>% group_by(Orthogroup, accession) %>% filter(all(selected_pfam_domains %in% PfamList)) %>% ungroup()
  ogs_long_select_pfam_has_domains_nodups <- ogs_long_select_pfam_has_domains %>% filter(!duplicated(accession))
  
  resampled_OGs_long_curr <- data.frame(accession = ogs_long_select_pfam_has_domains_nodups$accession, Orthogroup = OG_domains_curr$new_OG_id, taxid = ogs_long_select_pfam_has_domains_nodups$taxid, BOOL_PRIMARY_OG = ogs_long_select_pfam_has_domains_nodups$BOOL_PRIMARY_OG)
  resampled_OGs_long <- rbind(resampled_OGs_long, resampled_OGs_long_curr)
}

# Write out
write.table(resampled_OGs_long, here("data", "orthogroups", "refined_orthogroups", "refined_subsampled_OGs_euk203spp_long.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

