# Load libraries
library(here)
library(tidyverse)

## Read in datasets
# Read in DeepLoc results
deeploc_results <- read.table(here("data", "deeploc", "predictions", "DeepLoc2.0-mito_predictions.tsv"), header=TRUE)
colnames(deeploc_results) <- c("Protein_ID", "Mitochondrion")
deeploc_results$taxid <- gsub("_.*", "", deeploc_results$Protein_ID)

uniprot_proteomes_tax <- read.table(here("data", "taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
deeploc_results <- deeploc_results %>% filter(taxid %in% uniprot_proteomes_tax$TaxId)

deeploc_thresholds <- read.csv(here("data", "deeploc", "deeploc_thresholds.csv"))
deeploc_mito_threshold <- deeploc_thresholds$threshold[deeploc_thresholds$label == "Mitochondrion"]

deeploc_results$label <- "Nonmito"
deeploc_results$label[deeploc_results$Mitochondrion >= deeploc_mito_threshold] <- "Mito"

# Assign OGs
ogs_long <- read.table(here("data", "orthogroups", "refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")
ogs_long_primary <- ogs_long %>% filter(BOOL_PRIMARY_OG)
ogs_long_reduce <- ogs_long[,c("accession", "Orthogroup", "BOOL_PRIMARY_OG")]
deeploc_results <- merge(deeploc_results, ogs_long_reduce, by.x="Protein_ID", by.y="accession") # Allow non-primary OGs (many-to-many)

# Mark organelle-encoded proteins, if available
all_nonmito_accessions <- read.table(here("data", "deeploc", "all_nonmito_organelle_dna_protein_accessions_combined.txt"))$V1
all_mtdna_accessions <- read.table(here("data", "deeploc", "all_mtdna_protein_accessions_combined.txt"))$V1
deeploc_results$label[deeploc_results$Protein_ID %in% all_nonmito_accessions] <- "Nonmito"
deeploc_results$label[deeploc_results$Protein_ID %in% all_mtdna_accessions] <- "Mito"

# Ignore species missing organelle genomes in OGs that have organelle encoded proteins
missing_nonmito_organelle_taxids <- read.table(here("data", "deeploc", "missing_nonmito_organelle_taxids.txt"))$V1
missing_mtdna_taxids <- read.table(here("data", "deeploc", "missing_mtdna_taxids.txt"))$V1
nonmito_organelle_OG_ids <- unique(deeploc_results$Orthogroup[deeploc_results$Protein_ID %in% all_nonmito_accessions])
deeploc_results <- deeploc_results %>% filter(!(Orthogroup %in% nonmito_organelle_OG_ids & taxid %in% missing_nonmito_organelle_taxids))
mito_organelle_OG_ids <- unique(deeploc_results$Orthogroup[deeploc_results$Protein_ID %in% all_mtdna_accessions])
deeploc_results <- deeploc_results %>% filter(!(Orthogroup %in% mito_organelle_OG_ids & taxid %in% missing_mtdna_taxids))

deeploc_results$species <- uniprot_proteomes_tax$ScientificName[match(deeploc_results$taxid, uniprot_proteomes_tax$tree_id)]

gold_gene_accession_OG_id_df <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_2026.04.05.tsv"), sep="\t", header=TRUE)

deeploc_results$type <- "Predicted"
experimental_mitoproteome_species <- c("9606", "559292", "3702", "1257118", "5741", "5689", "185431", "32595")
deeploc_results$type[deeploc_results$taxid %in% experimental_mitoproteome_species] <- "Experimental"
deeploc_results$label[deeploc_results$taxid %in% experimental_mitoproteome_species] <- "Nonmito"
deeploc_results$label[deeploc_results$Protein_ID %in% gold_gene_accession_OG_id_df$gene_accession[gold_gene_accession_OG_id_df$taxid %in% experimental_mitoproteome_species]] <- "Mito"

gold_gene_accession_OG_id_df_primary <- read.table(here("data", "mito_orthogroups", "mito_proteins_experimental.and.mtDNA_primary.OG_2026.04.05.tsv"), sep="\t", header=TRUE)
primary_mito_OGs_mitocarta_or_mtdna <- unique(unlist(strsplit(gold_gene_accession_OG_id_df_primary$OG_id, split=",")))


## Process species mitoproteomes using species conservation filter
n_mito_species_threshold <- 5
deeploc_results$label <- factor(deeploc_results$label, levels = c("Nonmito", "Possible mito", "Mito", "MitoCarta", "Mito in <5 species"))

# Find OGs that have few species with a predicted mito protein. Use primary OGs to avoid double-counting proteins in each species
deeploc_results_lt5spp <- deeploc_results %>% filter(!Orthogroup %in% primary_mito_OGs_mitocarta_or_mtdna) %>% group_by(Orthogroup) %>% filter(label %in% c("Mito", "MitoCarta")) %>% filter(length(unique(taxid)) < n_mito_species_threshold)
deeploc_results_lt5spp_OG_ids <- unique(deeploc_results_lt5spp$Orthogroup)

# Label predictions using species conservation filter
deeploc_results_ogs_primary <- deeploc_results %>% filter(BOOL_PRIMARY_OG)
deeploc_results_ogs_primary$Orthogroup <- ogs_long_primary$Orthogroup[match(deeploc_results_ogs_primary$Protein_ID, ogs_long_primary$accession)]
deeploc_results_ogs_primary$label[deeploc_results_ogs_primary$Orthogroup %in% deeploc_results_lt5spp_OG_ids & deeploc_results_ogs_primary$label == "Mito"] <- "Mito in <5 species"

deeploc_results_ogs_primary_reduce <- deeploc_results_ogs_primary[,c("Protein_ID", "Orthogroup", "species", "taxid", "type", "label", "Mitochondrion")]
colnames(deeploc_results_ogs_primary_reduce) <- c("Protein_ID", "Orthogroup", "Species", "TaxID", "Proteome_type", "Prediction", "Predicted_mito_score")
write.table(deeploc_results_ogs_primary_reduce, here("deeploc_results_primary_OGs_processed.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

