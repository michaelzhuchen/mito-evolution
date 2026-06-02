---
  title: "R Notebook"
output: html_notebook
---
  
  This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 

```{r}
library(tidyverse)
```


```{r}
# ogs_long <- read.table("/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/uniprot_proteomes/OG_accessions/euk673spp_prokgroups_orthofinder_OG_all_merged.hmm.foldseek_multidomain_add.species.singleton.mtDNA_OG.UOG_protein_accessions_2025.05.31_long.txt", sep="\t", header=TRUE)
# ogs_long <- read.table("/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/uniprot_proteomes/OG_accessions/euk673spp_prokgroups_orthofinder_OG_all_merged.hmm.foldseek_multidomain_add.species.singleton.mtDNA_OG.UOG_protein_accessions_2025.08.12_long.txt", sep="\t", header=TRUE)
# ogs_long <- read.table("/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/uniprot_proteomes/OG_accessions/euk203spp_prokgroups_orthofinder_OG_all_merged.hmm.foldseek_multidomain_add.species.singleton.mtDNA_OG.UOG_protein_accessions_2025.05.31_long.txt", sep="\t", header=TRUE)
ogs_long <- read.table("/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/uniprot_proteomes/OG_accessions/euk203spp_prokgroups_orthofinder_OG_all_merged.hmm.foldseek_multidomain_add.species.singleton.mtDNA_OG.UOG_protein_accessions_2025.08.12_long.txt", sep="\t", header=TRUE)

# Get just mito OGs
# gold_gene_accession_OG_id_df <- read.table('/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/mito_lists/orthogroups/goldp_gene_to_OG.UOG_euk673spp_prokgroups_OG_all_merged.hmm.foldseek_multidomain_add.species.singleton.mtDNA_nomissing_2025.05.04.tsv', sep="\t", header=TRUE)
# ogs_long <- ogs_long %>% filter(Orthogroup %in% gold_gene_accession_OG_id_df$OG_id)

# To reduce computation, only examine OGs containing any MITOEPI species
eukaryote_reference_species_list <- c("9606", "559292", "3702", "1257118", "5689", "185431", "5741", "36329", "32595", "508771")
selected_OG_ids <- unique(ogs_long$Orthogroup[ogs_long$taxid %in% eukaryote_reference_species_list])
ogs_long <- ogs_long %>% filter(Orthogroup %in% selected_OG_ids)

uniprot_proteomes_tax <- readRDS("/Users/cmichael/Dropbox (MIT)/Mootha_lab/phylogenetics/data/uniprot_proteomes/uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.rds")
uniprot_proteomes_all_tax <- readRDS("/Users/cmichael/Dropbox (MIT)/Mootha_lab/phylogenetics/data/uniprot_proteomes/uniprot_bacteria.downsample.level6_tax_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.rds")

# Read in the prokaryote mmseqs2 cluster data
prokaryote_mmseqs2_clusters_taxids <- read.table('/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/uniprot_proteomes/prokgroups_mmseqs2_clusters/prokaryote_mmseqs2_clusters_taxids.tsv', sep="\t", quote="", header=TRUE)
```

```{r}
# Single species
ogs_long_single_species <- ogs_long %>% group_by(Orthogroup) %>% summarize(n_species = length(unique(taxid))) %>% filter(n_species == 1)
```


```{r}
# Prokaryotic

# Just use 2+ prok clades
ogs_long_prok <- ogs_long %>% filter(taxid %in% uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain != "Eukaryota"])
ogs_long_prokgroups <- ogs_long_prok %>% group_by(Orthogroup) %>% filter(length(unique(taxid)) > 1)

# Get prokaryote proteins and species counts
ogs_long_prok$prok_protein_id <- sub("^[^_]*_", "", ogs_long_prok$accession)
ogs_long_prok$prok_taxids <- prokaryote_mmseqs2_clusters_taxids$unique_member_seq_taxids[match(ogs_long_prok$prok_protein_id, prokaryote_mmseqs2_clusters_taxids$rep_seq)]

ogs_long_prok <- ogs_long_prok %>% separate_rows(prok_taxids, sep=",")
ogs_long_prok_count <- ogs_long_prok %>% group_by(Orthogroup) %>% summarize(n_prok_species = length(unique(prok_taxids)))

# # Use same parameters as infer origin for trees
# min_frac_prok_species_per_subgroup_for_prok_origin <- 0.25
# min_frac_subgroups_per_prokgroup_for_prok_origin <- 0.5
# min_prok_species_in_clade_to_not_prune <- 10
# 
# # Get total species counts per prokaryotic subgroup
# uniprot_proteomes_bacteria_subgroup_counts <- uniprot_proteomes_all_tax %>% filter(domain == "Bacteria") %>% group_by(tree_id, level_5) %>% summarize(n_species_per_subgroup = n(), .groups="keep") %>% mutate(n_species_per_subgroup_threshold = n_species_per_subgroup*min_frac_prok_species_per_subgroup_for_prok_origin)
# uniprot_proteomes_archaea_subgroup_counts <- uniprot_proteomes_all_tax %>% filter(domain == "Archaea") %>% group_by(tree_id, kingdom) %>% summarize(n_species_per_subgroup = n(), .groups="keep") %>% mutate(n_species_per_subgroup_threshold = n_species_per_subgroup*min_frac_prok_species_per_subgroup_for_prok_origin)
# colnames(uniprot_proteomes_bacteria_subgroup_counts)[2] <- "subgroup"
# colnames(uniprot_proteomes_archaea_subgroup_counts)[2] <- "subgroup"
# uniprot_proteomes_prokaryote_subgroup_counts <- rbind(uniprot_proteomes_bacteria_subgroup_counts, uniprot_proteomes_archaea_subgroup_counts)
# uniprot_proteomes_prokaryote_subgroup_counts_summary <- uniprot_proteomes_prokaryote_subgroup_counts %>% group_by(tree_id) %>% summarize(n_subgroups_per_prokgroup = n(), .groups = "keep") %>% mutate(n_subgroups_per_prokgroup_threshold = n_subgroups_per_prokgroup*min_frac_subgroups_per_prokgroup_for_prok_origin)
# 
# # Get species counts per prokaryotic subgroup for largest prokaryotic clade
# count_large_prokgroups_in_clade <- function(prok_species_largest_prok_clade) {
#   # prok_species_largest_prok_clade <- unique(unlist(strsplit(prok_species_in_clade, split=",")))
#   prok_largest_clade_tax <- uniprot_proteomes_all_tax[which(uniprot_proteomes_all_tax$TaxId %in% prok_species_largest_prok_clade),]
#   prok_largest_clade_bacteria_subgroup_counts <- prok_largest_clade_tax %>% filter(domain == "Bacteria") %>% group_by(tree_id, level_5) %>% summarize(n_species_per_subgroup = n(), .groups="keep")
#   prok_largest_clade_archaea_subgroup_counts <- prok_largest_clade_tax %>% filter(domain == "Archaea") %>% group_by(tree_id, kingdom) %>% summarize(n_species_per_subgroup = n(), .groups="keep")
#   colnames(prok_largest_clade_bacteria_subgroup_counts)[2] <- "subgroup"
#   colnames(prok_largest_clade_archaea_subgroup_counts)[2] <- "subgroup"
#   prok_largest_clade_prokaryote_subgroup_counts <- rbind(prok_largest_clade_bacteria_subgroup_counts, prok_largest_clade_archaea_subgroup_counts)
#   
#   # Identify prokaryotic subgroups passing the species thresholds
#   prok_largest_clade_prokaryote_subgroup_counts$n_species_per_subgroup_threshold <- uniprot_proteomes_prokaryote_subgroup_counts$n_species_per_subgroup_threshold[match(prok_largest_clade_prokaryote_subgroup_counts$subgroup, uniprot_proteomes_prokaryote_subgroup_counts$subgroup)]
#   prok_largest_clade_prokaryote_subgroup_counts_filter <- prok_largest_clade_prokaryote_subgroup_counts[which(prok_largest_clade_prokaryote_subgroup_counts$n_species_per_subgroup >= prok_largest_clade_prokaryote_subgroup_counts$n_species_per_subgroup_threshold),]
#   prok_largest_clade_prokaryote_prokgroup_subgroups <- prok_largest_clade_prokaryote_subgroup_counts_filter %>% group_by(tree_id) %>% summarize(n_subgroups_per_prokgroup=n())
#   
#   # Identify prokaryotic groups passing the subgroup thresholds
#   prok_largest_clade_prokaryote_prokgroup_subgroups$n_subgroups_per_prokgroup_threshold <- uniprot_proteomes_prokaryote_subgroup_counts_summary$n_subgroups_per_prokgroup_threshold[match(prok_largest_clade_prokaryote_prokgroup_subgroups$tree_id, uniprot_proteomes_prokaryote_subgroup_counts_summary$tree_id)]
#   prok_largest_clade_prokaryote_prokgroup_filter <- prok_largest_clade_prokaryote_prokgroup_subgroups[which(prok_largest_clade_prokaryote_prokgroup_subgroups$n_subgroups_per_prokgroup >= prok_largest_clade_prokaryote_prokgroup_subgroups$n_subgroups_per_prokgroup_threshold),]
#   
#   return(nrow(prok_largest_clade_prokaryote_prokgroup_filter))
#   # return(paste0(prok_largest_clade_prokaryote_prokgroup_filter$tree_id, collapse=","))
# }
# 
# ogs_long_prokgroups <- ogs_long_prok %>% group_by(Orthogroup) %>% filter(length(unique(prok_taxids)) >= min_prok_species_in_clade_to_not_prune) %>% summarize(n_large_prokgroups_in_clade = count_large_prokgroups_in_clade(prok_taxids))

# ogs_long_prokgroups <- ogs_long_prok %>% group_by(Orthogroup) %>% filter(length(unique(prok_taxids)) >= min_prok_species_in_clade_to_not_prune) %>% summarize(n_large_prokgroups_in_clade = count_large_prokgroups_in_clade(prok_taxids))
# ogs_long_prokgroups_seprows <- ogs_long_prokgroups %>% separate_rows(n_large_prokgroups_in_clade, sep=",")
# ogs_long_prokgroups_seprows <- ogs_long_prokgroups_seprows %>% filter(n_large_prokgroups_in_clade %in% large_bacterial_kingdoms)
# OG_ids_to_update <- unique(ogs_long_prokgroups_seprows$Orthogroup)
# write.table(OG_ids_to_update, "~/Downloads/OG_ids_to_update_infer_origin_prune_root_v6.4.1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

```{r}
## Eukaryotic
# Get eukaryote proteins and species counts. For eukaryotic, require at least 2 species in each of at least 2 eukaryotic clades
min_n_species_per_group <- 2

ogs_long_euk <- ogs_long %>% filter(taxid %in% uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Eukaryota"])
ogs_long_euk_count <- ogs_long_euk %>% group_by(Orthogroup) %>% summarize(n_euk_species = length(unique(taxid)))

ogs_long_species_count <- merge(ogs_long_euk_count, ogs_long_prok_count, by = "Orthogroup", all = TRUE)
ogs_long_species_count$n_euk_species[is.na(ogs_long_species_count$n_euk_species)] <- 0
ogs_long_species_count$n_prok_species[is.na(ogs_long_species_count$n_prok_species)] <- 0
ogs_long_species_count <- ogs_long_species_count %>% rowwise() %>% mutate(fraction_euk_species = n_euk_species / (n_euk_species + n_prok_species))
ogs_long_species_count <- ogs_long_species_count %>% filter(fraction_euk_species > 0.9)

ogs_long_eukgroups <- ogs_long_euk %>% filter(Orthogroup %in% ogs_long_species_count$Orthogroup)
ogs_long_eukgroups$euk_group <- uniprot_proteomes_tax$superfamily[match(ogs_long_eukgroups$taxid, uniprot_proteomes_tax$TaxId)]
ogs_long_eukgroups <- ogs_long_eukgroups %>% group_by(Orthogroup, euk_group) %>% mutate(n_species_per_group = n()) %>% filter(n_species_per_group >= min_n_species_per_group) %>% group_by(Orthogroup) %>% summarize(n_euk_groups = length(unique(euk_group)), euk_groups = paste0(unique(euk_group), collapse=","))

```


```{r}
# Assign origin
origin_df <- data.frame(Orthogroup = unique(ogs_long$Orthogroup), origin_label = NA)

single_species_OG_ids <- ogs_long_single_species$Orthogroup
origin_df$origin_label[origin_df$Orthogroup %in% single_species_OG_ids] <- "Species"

ogs_long_eukgroups_clade <- ogs_long_eukgroups[which(ogs_long_eukgroups$n_euk_groups == 1),]
origin_df$origin_label[match(ogs_long_eukgroups_clade$Orthogroup, origin_df$Orthogroup)] <- paste0("Clade_", ogs_long_eukgroups_clade$euk_groups)

eukaryota_OG_ids <- ogs_long_eukgroups$Orthogroup[which(ogs_long_eukgroups$n_euk_groups > 1)]
origin_df$origin_label[origin_df$Orthogroup %in% eukaryota_OG_ids] <- "Eukaryote"

# prokaryote_OG_ids <- ogs_long_prokgroups$Orthogroup[which(ogs_long_prokgroups$n_large_prokgroups_in_clade > 0)]
prokaryote_OG_ids <- ogs_long_prokgroups$Orthogroup
origin_df$origin_label[origin_df$Orthogroup %in% prokaryote_OG_ids] <- "Prokaryote"

origin_df$origin_label[which(is.na(origin_df$origin_label))] <- "Indeterminate"

sort(table(origin_df$origin_label), decreasing=TRUE)

write.table(origin_df, '/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/mitosynthesis/origin/OG_origin/euk203spp_prokgroups_OG_origin.tsv', sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
# write.table(origin_df, '/Users/cmichael/MIT Dropbox/Michael Chen/Mootha_lab/phylogenetics/data/mitosynthesis/origin/OG_origin/euk673spp_prokgroups_OG_origin.tsv', sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
```



Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

