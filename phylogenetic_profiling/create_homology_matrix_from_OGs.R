### Generate reference-free multispecies CLIME matrix directly from OGs for 129 eukaryotes that have a complete mtDNA annotation

# Load libraries
library(tidyverse)

# Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

core_species <- c("9606", "559292", "1257118", "5741", "5689", "185431", "3702", "32595")
eukaryote_reference_species_list <- c(core_species, uniprot_proteomes_tax$TaxId[uniprot_proteomes_tax$domain == "Eukaryota"])
eukaryote_reference_species_list <- eukaryote_reference_species_list[!duplicated(eukaryote_reference_species_list)]

prokaryote_tree_ids <- uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain != "Eukaryota"]

ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")

ogs_long_select_only <- ogs_long %>% filter(taxid %in% eukaryote_reference_species_list)


ogs_long_select_only_OGs <- unique(ogs_long_select_only$Orthogroup)
ogs_long_select <- ogs_long %>% filter(Orthogroup %in% ogs_long_select_only_OGs)

## Add collapsed prokaryote outgroup
# Get list of prokaryotic groups
ogs_long_outgroup <- ogs_long_select %>% filter(taxid %in% prokaryote_tree_ids) %>% group_by(Orthogroup) %>% filter(!duplicated(taxid))
# Filter for 2+ prokaryotic groups
ogs_long_outgroup <- ogs_long_outgroup %>% group_by(Orthogroup) %>% filter(length(unique(taxid)) > 1)
ogs_long_outgroup$taxid <- "Asgardgroup"
ogs_long_outgroup <- ogs_long_outgroup %>% group_by(Orthogroup) %>% filter(!duplicated(taxid))
# Remove and re-add prokaryote
ogs_long_select <- ogs_long_select %>% filter(!taxid %in% prokaryote_tree_ids)
ogs_long_select <- rbind(ogs_long_select, ogs_long_outgroup)

ogs_long_select <- ogs_long_select[,c("Orthogroup", "taxid")]
ogs_long_select_mat <- ogs_long_select %>%
  mutate(present = 1) %>%  # Add a binary column indicating presence
  distinct() %>%  # Remove duplicates, if any
  pivot_wider(
    names_from = taxid,
    values_from = present,
    values_fill = list(present = 0)  # Fill missing values with 0
  )

# Only use species that have a complete mtDNA annotation
missing_mtdna_taxids <- read.table(here("data/deeploc", "missing_mtdna_taxids.txt"))$V1
naegleria_lovaniensis_taxid <- 51637 # Exclude Naegleria lovaniensis due to incomplete mtDNA annotation
missing_mtdna_taxids <- c(missing_mtdna_taxids, naegleria_lovaniensis_taxid) 
ogs_long_select_mat <- ogs_long_select_mat[,which(!colnames(ogs_long_select_mat) %in% missing_mtdna_taxids)]

# Get human symbol mapping
# Read in mapping of human gene IDs to uniprot
human_id2training <- read.table(here("data/orthogroups/idmapping", "human.id2training_uniprot_symbols.txt"), sep="\t", header=TRUE)
human_id2training$accession <- paste0("9606_", human_id2training$Entry)
ogs_long_select_only$Entrez <- human_id2training$HumanGeneID[match(ogs_long_select_only$accession, human_id2training$accession)]

ogs_long_select_only$Symbol <- human_id2training$Symbol[match(ogs_long_select_only$accession, human_id2training$accession)]

# Get human gene symbols from CLIME matrix, if available
clime_mat <- read.table(here("data/phylogenetic_profiling/BLAST_euk138spp_prokaryote", "hsa.matrix138.e3.q00.p20.txt"), sep="\t", header=TRUE)
clime_mat_nomissingsymbol <- clime_mat[which(clime_mat$Symbol != ""),]
ogs_long_select_only$Symbol[ogs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez] <- clime_mat_nomissingsymbol$Symbol[match(ogs_long_select_only$Entrez[ogs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez], clime_mat_nomissingsymbol$Entrez)]
ogs_long_select_only$accession[ogs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez] <- paste0("hsa_", clime_mat_nomissingsymbol$Entrez[match(ogs_long_select_only$Entrez[ogs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez], clime_mat_nomissingsymbol$Entrez)])

# Get symbols for non-human species
id_mapping <- read.delim(here("data/orthogroups/idmapping", "mitoepi_10spp_geneid_to_symbol_no_version.tsv"), header=FALSE)
colnames(id_mapping) <- c("gene_id", "symbol", "taxid")
id_mapping <- id_mapping %>% rowwise() %>% mutate(accession = paste0(taxid, "_", gene_id))
id_mapping <- id_mapping %>% filter(!is.na(symbol))
ogs_long_select_only$Symbol[ogs_long_select_only$accession %in% id_mapping$accession] <- id_mapping$symbol[match(ogs_long_select_only$accession[ogs_long_select_only$accession %in% id_mapping$accession], id_mapping$accession)]

# Use accessions
taxids_order <- data.frame(taxid = eukaryote_reference_species_list)
taxids_order$index <- 1:nrow(taxids_order)
ogs_long_select_only$order_index <- taxids_order$index[match(ogs_long_select_only$taxid, taxids_order$taxid)]
ogs_long_select_only <- ogs_long_select_only[order(ogs_long_select_only$Orthogroup, is.na(ogs_long_select_only$Symbol), ogs_long_select_only$order_index, ogs_long_select_only$Symbol),]
ogs_long_select_only$Symbol[is.na(ogs_long_select_only$Symbol)] <- sub("^[^_]*_", "", ogs_long_select_only$accession[is.na(ogs_long_select_only$Symbol)])
ogs_long_select_only$Symbol <- gsub("-t26.*", "", ogs_long_select_only$Symbol)
ogs_long_select_only$Symbol <- gsub("PeptideAtlas_*", "", ogs_long_select_only$Symbol)
ogs_long_select_only$Symbol <- gsub("\\.1-p.*", "", ogs_long_select_only$Symbol)

ogs_long_select_only_symbol <- ogs_long_select_only %>% group_by(Orthogroup) %>% summarize(Symbol = paste0(Symbol, collapse=","))
ogs_long_select_only_symbol <- ogs_long_select_only_symbol %>% rowwise() %>% mutate(Symbol = gsub("(NA)", "", Symbol, fixed=TRUE))
ogs_long_select_only_symbol <- ogs_long_select_only_symbol %>% rowwise() %>% mutate(Symbol = paste0(Symbol, ";", Orthogroup))
n_max_characters <- 19
ogs_long_select_only_symbol <- ogs_long_select_only_symbol %>% rowwise() %>% mutate(Symbol = str_trunc(Symbol, n_max_characters, ellipsis="."))
ogs_long_select_only_symbol$Symbol <- make.unique(ogs_long_select_only_symbol$Symbol, sep="")

ogs_long_select_mat$Symbol <- ogs_long_select_only_symbol$Symbol[match(ogs_long_select_mat$Orthogroup, ogs_long_select_only_symbol$Orthogroup)]

ogs_long_select_mat$Entrez <- ogs_long_select_mat$Orthogroup

species_taxids <- sort(colnames(ogs_long_select_mat)[which(!colnames(ogs_long_select_mat) %in% c("Entrez", "Symbol", "Orthogroup"))])
ogs_long_select_mat_for_clime <- ogs_long_select_mat[,c("Entrez", "Symbol", species_taxids)]

# Rename to G.species
# Extract first letter of genus and species epithet
convert_to_abbreviation <- function(name) {
  # Split into words
  parts <- unlist(strsplit(name, "\\s+"))
  if (length(parts) < 2) return(name)  # Return name if not enough parts
  genus_abbr <- parts[1]
  species <- parts[2]
  strain <- parts[3]
  substrain <- parts[4]
  return(paste0(genus_abbr, ".", species, ".", strain, ".", substrain))
}
uniprot_proteomes_tax$species_name_for_clime <- sapply(uniprot_proteomes_tax$ScientificName, convert_to_abbreviation)
uniprot_proteomes_tax$species_name_for_clime <- gsub(".NA", "", uniprot_proteomes_tax$species_name_for_clime, fixed=TRUE)
uniprot_proteomes_tax$species_name_for_clime <- gsub("[^a-zA-Z0-9]", ".", uniprot_proteomes_tax$species_name_for_clime)
# Sanity check: no species should be duplicated
uniprot_proteomes_tax$species_name_for_clime[which(duplicated(uniprot_proteomes_tax$species_name_for_clime))]

# Rename columns to species names
colnames(ogs_long_select_mat_for_clime)[colnames(ogs_long_select_mat_for_clime) %in% uniprot_proteomes_tax$tree_id] <- uniprot_proteomes_tax$species_name_for_clime[match(colnames(ogs_long_select_mat_for_clime)[colnames(ogs_long_select_mat_for_clime) %in% uniprot_proteomes_tax$tree_id], uniprot_proteomes_tax$tree_id)]
# Rename prokaryotic outgroup
colnames(ogs_long_select_mat_for_clime)[which(colnames(ogs_long_select_mat_for_clime) == "Asgard.group")] <- "Prokaryotes"

## Write out
write.table(ogs_long_select_mat_for_clime, here("data/phylogenetic_profiling/OG_euk129spp.mtDNA_prokaryote", "OG.UOG_euk129spp_prokaryote_all_euks_mat_accessions_parent.OG.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

