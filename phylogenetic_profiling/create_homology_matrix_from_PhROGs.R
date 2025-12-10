### Generate reference-free multispecies CLIME matrix from PhROGs for 203 eukaryotes

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
progs_long <- read.table(file.path("data/phylogenetically_resolved_orthogroups", "PhROGs_long", paste0("PhROGs_at_Node34_Eukaryota_long.tsv")), sep="\t", header=TRUE)
colnames(progs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
progs_long <- progs_long %>% mutate(Orthogroup = PROG_id, OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id))
progs_long <- progs_long %>% group_by(OG_id) %>% filter(!duplicated(protein_id)) # assign each protein to its largest PhROG

progs_long_select_only <- progs_long %>% filter(taxid %in% eukaryote_reference_species_list)
progs_long_select_only_PROGs <- unique(progs_long_select_only$PROG_id)
progs_long_select_only_OGs <- unique(gsub("_.*", "", progs_long_select_only_PROGs))
progs_long_select <- progs_long %>% filter(PROG_id %in% progs_long_select_only_PROGs)
progs_long_select <- progs_long_select[,c("Orthogroup", "taxid")]

# Add collapsed prokaryotic outgroup
progs_long_select_for_merge <- progs_long_select %>% filter(!duplicated(Orthogroup))
progs_long_select_for_merge$OG_id <- gsub("_.*", "", progs_long_select_for_merge$Orthogroup)
progs_long_select_for_merge <- progs_long_select_for_merge[,c("Orthogroup", "OG_id")]
# Get list of prokaryotic groups
ogs_long_outgroup <- ogs_long %>% filter(taxid %in% prokaryote_tree_ids) %>% filter(Orthogroup %in% progs_long_select_for_merge$OG_id) %>% group_by(Orthogroup) %>% filter(!duplicated(taxid))

# Filter for 2+ prokaryotic groups
ogs_long_outgroup <- ogs_long_outgroup %>% group_by(Orthogroup) %>% filter(length(unique(taxid)) > 1)
ogs_long_outgroup$taxid <- "Asgardgroup"

ogs_long_outgroup <- ogs_long_outgroup %>% group_by(Orthogroup) %>% filter(!duplicated(taxid))
ogs_long_outgroup_merge <- merge(progs_long_select_for_merge, ogs_long_outgroup, by.x="OG_id", by.y="Orthogroup")
ogs_long_outgroup_merge <- ogs_long_outgroup_merge[,c("Orthogroup", "taxid")]
progs_long_select <- rbind(progs_long_select, ogs_long_outgroup_merge)


progs_long_select_mat <- progs_long_select %>%
  mutate(present = 1) %>%  # Add a binary column indicating presence
  distinct() %>%  # Remove duplicates, if any
  pivot_wider(
    names_from = taxid,
    values_from = present,
    values_fill = list(present = 0)  # Fill missing values with 0
  )


# Include only reference species OGs
ogs_long_select_only <- ogs_long %>% filter(!Orthogroup %in% progs_long_select_only_OGs) %>% filter(taxid %in% eukaryote_reference_species_list)

ogs_long_select_only_OGs <- unique(ogs_long_select_only$Orthogroup)
ogs_long_select <- ogs_long %>% filter(Orthogroup %in% ogs_long_select_only_OGs)

# Add collapsed prokaryotic outgroup
# Get list of prokaryotic groups
ogs_long_outgroup <- ogs_long_select %>% filter(taxid %in% prokaryote_tree_ids) %>% group_by(Orthogroup) %>% filter(!duplicated(taxid))

# Filter for 2+ prokaryotic groups
ogs_long_outgroup <- ogs_long_outgroup %>% group_by(Orthogroup) %>% filter(length(unique(taxid)) > 1)
ogs_long_outgroup$taxid <- "Asgardgroup"

ogs_long_outgroup <- ogs_long_outgroup %>% group_by(Orthogroup) %>% filter(!duplicated(taxid))
# Remove and re-add prokaryote
ogs_long_select <- ogs_long_select %>% filter(!taxid %in% prokaryote_tree_ids)
ogs_long_select <- rbind(ogs_long_select, ogs_long_outgroup)

ogs_long_select <- ogs_long_select %>% filter(taxid %in% colnames(progs_long_select_mat))
ogs_long_select <- ogs_long_select[,c("Orthogroup", "taxid")]
ogs_long_select_mat <- ogs_long_select %>%
  mutate(present = 1) %>%  # Add a binary column indicating presence
  distinct() %>%  # Remove duplicates, if any
  pivot_wider(
    names_from = taxid,
    values_from = present,
    values_fill = list(present = 0)  # Fill missing values with 0
  )
# Add any missing species columns
missing_columns_to_add <- colnames(progs_long_select_mat)[!colnames(progs_long_select_mat) %in% colnames(ogs_long_select_mat)]
ogs_long_select_mat_missing_columns_to_add <- matrix(data = 0, nrow = nrow(ogs_long_select_mat), ncol = length(missing_columns_to_add))
colnames(ogs_long_select_mat_missing_columns_to_add) <- missing_columns_to_add
ogs_long_select_mat <- cbind(ogs_long_select_mat, ogs_long_select_mat_missing_columns_to_add)

progs_long_select_mat <- rbind(progs_long_select_mat, ogs_long_select_mat)
progs_long_select_only <- progs_long_select_only[,c("Orthogroup", "protein_id", "taxid")]
ogs_long_select_only$protein_id <- ogs_long_select_only$accession
ogs_long_select_only <- ogs_long_select_only[,c("Orthogroup", "protein_id", "taxid")]
progs_long_select_only <- rbind(progs_long_select_only, ogs_long_select_only)


# Get human symbol mapping
human_id2training <- read.table(here("data/orthogroups/idmapping", "human.id2training_uniprot_symbols.txt"), sep="\t", header=TRUE)
human_id2training$accession <- paste0("9606_", human_id2training$Entry)
progs_long_select_only$Entrez <- human_id2training$HumanGeneID[match(progs_long_select_only$protein_id, human_id2training$accession)]

progs_long_select_only$Symbol <- human_id2training$Symbol[match(progs_long_select_only$protein_id, human_id2training$accession)]

# Get human gene symbols from CLIME matrix, if available
clime_mat <- read.table(here("data/phylogenetic_profiling/BLAST_euk138spp_prokaryote", "hsa.matrix138.e3.q00.p20.txt"), sep="\t", header=TRUE)
clime_mat_nomissingsymbol <- clime_mat[which(clime_mat$Symbol != ""),]
progs_long_select_only$Symbol[progs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez] <- clime_mat_nomissingsymbol$Symbol[match(progs_long_select_only$Entrez[progs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez], clime_mat_nomissingsymbol$Entrez)]
progs_long_select_only$protein_id[progs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez] <- paste0("hsa_", clime_mat_nomissingsymbol$Entrez[match(progs_long_select_only$Entrez[progs_long_select_only$Entrez %in% clime_mat_nomissingsymbol$Entrez], clime_mat_nomissingsymbol$Entrez)])


# Get symbols for non-human species
id_mapping <- read.delim(here("data/orthogroups/idmapping", "mitoepi_10spp_geneid_to_symbol_no_version.tsv"), header=FALSE)
colnames(id_mapping) <- c("gene_id", "symbol", "taxid")
id_mapping <- id_mapping %>% rowwise() %>% mutate(accession = paste0(taxid, "_", gene_id))
id_mapping <- id_mapping %>% filter(!is.na(symbol))
progs_long_select_only$Symbol[progs_long_select_only$protein_id %in% id_mapping$accession] <- id_mapping$symbol[match(progs_long_select_only$protein_id[progs_long_select_only$protein_id %in% id_mapping$accession], id_mapping$accession)]

# Use accessions
taxids_order <- data.frame(taxid = eukaryote_reference_species_list)
taxids_order$index <- 1:nrow(taxids_order)
progs_long_select_only$order_index <- taxids_order$index[match(progs_long_select_only$taxid, taxids_order$taxid)]
progs_long_select_only <- progs_long_select_only[order(progs_long_select_only$Orthogroup, is.na(progs_long_select_only$Symbol), progs_long_select_only$order_index, progs_long_select_only$Symbol),]
progs_long_select_only$Symbol[is.na(progs_long_select_only$Symbol)] <- sub("^[^_]*_", "", progs_long_select_only$protein_id[is.na(progs_long_select_only$Symbol)])
progs_long_select_only$Symbol <- gsub("-t26.*", "", progs_long_select_only$Symbol)
progs_long_select_only$Symbol <- gsub("PeptideAtlas_*", "", progs_long_select_only$Symbol)
progs_long_select_only$Symbol <- gsub("\\.1-p.*", "", progs_long_select_only$Symbol)

progs_long_select_only_symbol <- progs_long_select_only %>% group_by(Orthogroup) %>% summarize(Symbol = paste0(Symbol, collapse=","))
progs_long_select_only_symbol <- progs_long_select_only_symbol %>% rowwise() %>% mutate(Symbol = gsub("(NA)", "", Symbol, fixed=TRUE))
progs_long_select_only_symbol <- progs_long_select_only_symbol %>% rowwise() %>% mutate(Symbol = paste0(Symbol, ";", Orthogroup))
n_max_characters <- 19
progs_long_select_only_symbol <- progs_long_select_only_symbol %>% rowwise() %>% mutate(Symbol = str_trunc(Symbol, n_max_characters, ellipsis="."))
progs_long_select_only_symbol$Symbol <- make.unique(progs_long_select_only_symbol$Symbol, sep="")

progs_long_select_mat$Symbol <- progs_long_select_only_symbol$Symbol[match(progs_long_select_mat$Orthogroup, progs_long_select_only_symbol$Orthogroup)]

progs_long_select_mat$Entrez <- progs_long_select_mat$Orthogroup

species_taxids <- sort(colnames(progs_long_select_mat)[which(!colnames(progs_long_select_mat) %in% c("Entrez", "Symbol", "Orthogroup"))])
progs_long_select_mat_for_clime <- progs_long_select_mat[,c("Entrez", "Symbol", species_taxids)]

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
# Sanity check: none should be duplicated
uniprot_proteomes_tax$species_name_for_clime[which(duplicated(uniprot_proteomes_tax$species_name_for_clime))]

# Rename columns to species names
colnames(progs_long_select_mat_for_clime)[colnames(progs_long_select_mat_for_clime) %in% uniprot_proteomes_tax$tree_id] <- uniprot_proteomes_tax$species_name_for_clime[match(colnames(progs_long_select_mat_for_clime)[colnames(progs_long_select_mat_for_clime) %in% uniprot_proteomes_tax$tree_id], uniprot_proteomes_tax$tree_id)]
# Rename prokaryotic outgroup
colnames(progs_long_select_mat_for_clime)[which(colnames(progs_long_select_mat_for_clime) == "Asgard.group")] <- "Prokaryotes"

## Write out
write.table(progs_long_select_mat_for_clime, here("data/phylogenetic_profiling/PhROG_euk203spp_prokaryote", "OG.PhROG.Eukaryota_euk203spp_prokaryote_all_euks_mat_accessions.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
