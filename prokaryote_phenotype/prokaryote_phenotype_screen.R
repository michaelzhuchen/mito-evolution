### Screen for prokaryote phenotype enrichment 

# Load libraries
library(tidyverse)
library(here)

### Get prokaryote species phenotypes
uniprot_proteomes_all_tax <- read.table(here("data/taxonomy", "uniprot_bacteria.downsample.level6_tax_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

# read in GOLD species phenotypes
GOLD_phenotypes <- read.csv("~/Dropbox (MIT)/Mootha_lab/reference/GOLD/goldData_Organism.csv", header=TRUE)
GOLD_phenotypes_selected <- GOLD_phenotypes[which(GOLD_phenotypes$ORGANISM.NCBI.TAX.ID %in% uniprot_proteomes_all_tax$TaxId),]
GOLD_phenotypes_selected_nodups <- GOLD_phenotypes_selected[!duplicated(GOLD_phenotypes_selected$ORGANISM.NCBI.TAX.ID),]
uniprot_proteomes_all_tax_GOLD <- merge(uniprot_proteomes_all_tax, GOLD_phenotypes_selected_nodups, by.x="TaxId", by.y="ORGANISM.NCBI.TAX.ID", all.x=TRUE, all.y=FALSE)
length(which(uniprot_proteomes_all_tax_GOLD$ORGANISM.OXYGEN.REQUIREMENT != "" & uniprot_proteomes_all_tax_GOLD$domain == "Bacteria"))

# Import bacdive data
# read in the modified bacteria taxonomic info spreadsheet
bacdive_phenotypes <- read.csv(here("data/prokaryote_phenotype", "bacdive_phenotypes.csv"), header=TRUE)
# remove missing and duplicates
bacdive_phenotypes_tax_nomissing <- bacdive_phenotypes
bacdive_phenotypes_tax_nomissing <- bacdive_phenotypes_tax_nomissing[!duplicated(bacdive_phenotypes_tax_nomissing$TaxId),]
# manually correct incorrect entries
bacdive_phenotypes_tax_nomissing$oxygen_tolerance[bacdive_phenotypes_tax_nomissing$oxygen_tolerance == "25036"] <- "microaerophile"

# merge GOLD and bacdive phenotypes
uniprot_proteomes_all_tax_GOLD_bacdive <- merge(uniprot_proteomes_all_tax_GOLD, bacdive_phenotypes_tax_nomissing, by.x="TaxId", by.y="TaxId", all.x=TRUE, all.y=FALSE, no.dups=TRUE)
uniprot_proteomes_all_tax_GOLD_bacdive[duplicated(uniprot_proteomes_all_tax_GOLD_bacdive$TaxId),]

uniprot_proteomes_all_tax_GOLD_bacdive$ORGANISM.OXYGEN.REQUIREMENT <- tolower(uniprot_proteomes_all_tax_GOLD_bacdive$ORGANISM.OXYGEN.REQUIREMENT)

# determine which set of annotations for oxygen requirement are more trustworthy
uniprot_proteomes_all_tax_GOLD_bacdive_compare <- uniprot_proteomes_all_tax_GOLD_bacdive %>% filter(ORGANISM.OXYGEN.REQUIREMENT != "" & oxygen_tolerance != "") %>% filter(ORGANISM.OXYGEN.REQUIREMENT != oxygen_tolerance)
uniprot_proteomes_all_tax_GOLD_bacdive_compare[,c("ScientificName", "ORGANISM.OXYGEN.REQUIREMENT", "oxygen_tolerance")]

# Based on manual curation, if have two sets of conflicting annotations, go with the bacdive data in most cases, with a few exceptions that are manually corrected
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot <- uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_tolerance
missing_oxygen_annot <- which(is.na(uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot) | uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot == "")
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[missing_oxygen_annot] <- uniprot_proteomes_all_tax_GOLD_bacdive$ORGANISM.OXYGEN.REQUIREMENT[missing_oxygen_annot]
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive$ScientificName == "Thalassospira profundimaris"] <- "aerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive$ScientificName == "Haematobacter massiliensis"] <- "aerobe"

# rename categories
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot == "microaerophilic"] <- "microaerophile"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot == "obligate anaerobe"] <- "anaerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot == "obligate aerobe"] <- "aerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(is.na(uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot))] <- ""

# for eukaryotes, manually add in annotations
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$domain == "Eukaryota")] <- "aerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$level_6 == "Microsporidia")] <- "anaerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$superfamily == "Metamonada")] <- "anaerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$level_7 == "Myxozoa")] <- "anaerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$level_11 == "Cryptosporidium")] <- "anaerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$level_9 == "Neocallimastigales")] <- "anaerobe"
uniprot_proteomes_all_tax_GOLD_bacdive$oxygen_combined_annot[which(uniprot_proteomes_all_tax_GOLD_bacdive$ScientificName == "Entamoeba histolytica")] <- "anaerobe"


### Read in orthogroup data
ogs_long <- read.table(here("data/orthogroups", "refined_OGs_euk673spp_long.txt"), sep="\t", header=TRUE)

origin_table <- read.table(here("data/protein_phylogeny", "orthogroup_origin_domain.tsv"), sep="\t")
colnames(origin_table) <- c("OG_id", "n_euk_species_largest_euk_clade", "n_prok_species_largest_prok_clade", "LCA_node", "LCA_node_euks", "origin_domain", "dropped_tips")
origin_table$OG_id <- gsub(".faa_clipkit.gappy.msa", "", origin_table$OG_id, fixed=TRUE)
prok_origin_OG_ids <- origin_table$OG_id[origin_table$origin_domain == "Prokaryote"]

# Read in the prokaryote mmseqs2 cluster data
prokaryote_mmseqs2_clusters_taxids <- read.table(here("data/downsample_prokaryotes", "prokaryote_mmseqs2_clusters_taxids.tsv"), sep="\t", quote="", header=TRUE)

# Read in experimental and mtDNA mito proteins
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2025.09.30.tsv"), sep="\t", header=TRUE)

### Screen for enrichment
ogs_long_prok <- ogs_long %>% filter(Orthogroup %in% prok_origin_OG_ids) %>% filter(taxid %in% uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain != "Eukaryota"])

ogs_long_prok$prok_protein_id <- sub("^[^_]*_", "", ogs_long_prok$accession)

ogs_long_prok$prok_taxids <- prokaryote_mmseqs2_clusters_taxids$unique_member_seq_taxids[match(ogs_long_prok$prok_protein_id, prokaryote_mmseqs2_clusters_taxids$rep_seq)]

ogs_long_prok_summary <- ogs_long_prok %>% group_by(Orthogroup) %>% summarize(prok_taxids = paste0(unique(unlist(strsplit(prok_taxids, split=","))), collapse=","))

# Use GOLD + BacDive
uniprot_proteomes_all_tax_GOLD_bacdive_prok <- uniprot_proteomes_all_tax_GOLD_bacdive %>% filter(domain != "Eukaryota")
ogs_long_prok_summary_phenotypes <- ogs_long_prok_summary %>% rowwise() %>% mutate(n_aerobe = sum(uniprot_proteomes_all_tax_GOLD_bacdive_prok$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive_prok$TaxId %in% unlist(strsplit(prok_taxids, split=","))] %in% c("aerobe")), n_anaerobe = sum(uniprot_proteomes_all_tax_GOLD_bacdive_prok$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive_prok$TaxId %in% unlist(strsplit(prok_taxids, split=","))] %in% c("anaerobe", "microaerophile")), n_facultative = sum(uniprot_proteomes_all_tax_GOLD_bacdive_prok$oxygen_combined_annot[uniprot_proteomes_all_tax_GOLD_bacdive_prok$TaxId %in% unlist(strsplit(prok_taxids, split=","))] %in% c("facultative", "facultative aerobe", "facultative anaerobe")), n_gram_negative = sum(uniprot_proteomes_all_tax_GOLD_bacdive_prok$ORGANISM.GRAM.STAIN[uniprot_proteomes_all_tax_GOLD_bacdive_prok$TaxId %in% unlist(strsplit(prok_taxids, split=","))] %in% c("Gram-")), n_gram_positive = sum(uniprot_proteomes_all_tax_GOLD_bacdive_prok$ORGANISM.GRAM.STAIN[uniprot_proteomes_all_tax_GOLD_bacdive_prok$TaxId %in% unlist(strsplit(prok_taxids, split=","))] %in% c("Gram+")))
prokaryote_phenotypes_total <- uniprot_proteomes_all_tax_GOLD_bacdive_prok %>% filter(TaxId %in% unique(unlist(strsplit(ogs_long_prok_summary$prok_taxids, split=",")))) %>% summarize(n_aerobe_total = sum(oxygen_combined_annot %in% c("aerobe")), n_anaerobe_total = sum(oxygen_combined_annot %in% c("anaerobe", "microaerophile")), n_facultative_total = sum(oxygen_combined_annot %in% c("facultative", "facultative aerobe", "facultative anaerobe")), n_gram_negative_total = sum(ORGANISM.GRAM.STAIN %in% c("Gram-")), n_gram_positive_total = sum(ORGANISM.GRAM.STAIN %in% c("Gram+")))

ogs_long_prok_summary_phenotypes <- ogs_long_prok_summary_phenotypes %>% rowwise() %>% mutate(prop_aerobe = n_aerobe / prokaryote_phenotypes_total$n_aerobe_total, prop_anaerobe = n_anaerobe / prokaryote_phenotypes_total$n_anaerobe_total, prop_facultative = n_facultative / prokaryote_phenotypes_total$n_facultative_total, prop_gram_negative = n_gram_negative / prokaryote_phenotypes_total$n_gram_negative_total, prop_gram_positive = n_gram_positive / prokaryote_phenotypes_total$n_gram_positive_total)

# Filter for LECA mito
leca_PhROG_OG_ids <- read.table(here("data/reconstruction", 'leca_PhROG_OG_ids_corespecies_deeplocmc.min5spp_Eukaryota.parent.2supergroups.min4spp.filter.txt'))$V1
LECA_OG_ids <- unique(gsub("_.*", "", leca_PhROG_OG_ids))
leca_mito_PhROG_ids <- read.table(here("data/reconstruction", 'leca_mito_PhROG_ids_corespecies_deeplocmc.min5spp_Eukaryota.parent.2supergroups.min4spp.filter.txt'))$V1
LECA_mito_OG_ids <- unique(gsub("_.*", "", leca_mito_PhROG_ids))
ogs_long_prok_summary_phenotypes <- ogs_long_prok_summary_phenotypes %>% mutate(in_LECA = Orthogroup %in% LECA_OG_ids, in_LECA_mito = Orthogroup %in% LECA_mito_OG_ids)
# Select LECA +/- mito
ogs_long_prok_summary_phenotypes <- ogs_long_prok_summary_phenotypes %>% filter(in_LECA_mito)

# Compute log odds
ogs_long_prok_summary_phenotypes$log_odds_aerobe_vs_anaerobe <- log(ogs_long_prok_summary_phenotypes$prop_aerobe / ogs_long_prok_summary_phenotypes$prop_anaerobe)
ogs_long_prok_summary_phenotypes$log_odds_anaerobe_vs_aerobe <- log(ogs_long_prok_summary_phenotypes$prop_anaerobe / ogs_long_prok_summary_phenotypes$prop_aerobe)
ogs_long_prok_summary_phenotypes$log_odds_gramnegative_vs_grampositive <- log(ogs_long_prok_summary_phenotypes$prop_gram_negative / ogs_long_prok_summary_phenotypes$prop_gram_positive)

ogs_long_prok_summary_phenotypes <- ogs_long_prok_summary_phenotypes[order(ogs_long_prok_summary_phenotypes$log_odds_aerobe_vs_anaerobe, decreasing=FALSE),]

ogs_long_prok_summary_phenotypes$Rank <- 1:nrow(ogs_long_prok_summary_phenotypes)



