### HGT analysis

# Load libraries
library(tidyverse)
library(ape)
library(castor)
library(phytools)
library(RColorBrewer)

dataset_name <- "species_tree_1"

## Read in data
# Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data", "taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)
uniprot_proteomes_all_tax <- read.delim(here("data", "taxonomy", "uniprot_bacteria.downsample.level6_tax_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

# Read in Eukaryota parent PhROGs
parent_progs_long <- read.table(here("data", "phylogenetically_resolved_orthogroups", dataset_name, "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(parent_progs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")
parent_progs_long <- parent_progs_long %>% mutate(OG_id = gsub("_.*", "", PROG_id), taxid = gsub("_.*", "", protein_id)) %>% group_by(OG_id) %>% filter(!duplicated(protein_id))

# Get mito Eukaryota_parent PhROGs
parent_progs_long_mito_PhROG_ids <- read.table(here("data", "reconstruction", dataset_name, "parent_mito_PhROG_ids.txt"))$V1

# Read in HGT PhROGs
combined_hgt_phrog_raw <- read.table(here("data", "horizontal_gene_transfer", dataset_name, "posterior_clades_HGT_Node34_Eukaryota_parent_wide.tsv"), sep="\t")
colnames(combined_hgt_phrog_raw) <- c("OG_id","label","clade_index","distance_to_root","reference_protein_ids","count","species_overlap","duplications_rec","nonvertical_protein_ids","n_species","n_reference_proteins","fraction_euk_species","mito_localization_prob","PhROG_id","fraction_primary_OG_for_vertical_proteins","self_clade_protein_ids","sister_clade_protein_ids","cousin_clade_protein_ids","grandma_clade_protein_ids","self_clade_prok_species_in_clades","self_clade_n_prok_species","self_clade_n_euk_species","self_clade_fraction_euk_species","sister_clade_prok_species_in_clades","sister_clade_n_prok_species","sister_clade_n_euk_species","sister_clade_fraction_euk_species","cousin_clade_prok_species_in_clades","cousin_clade_n_prok_species","cousin_clade_n_euk_species","cousin_clade_fraction_euk_species","grandma_clade_prok_species_in_clades","grandma_clade_n_prok_species","grandma_clade_n_euk_species","grandma_clade_fraction_euk_species","HGT_self", "HGT_sister","HGT_cousin")

# Convert NAs to empty string to avoid errors
combined_hgt_phrog_raw$self_clade_protein_ids[is.na(combined_hgt_phrog_raw$self_clade_protein_ids)] <- ""
combined_hgt_phrog_raw$sister_clade_protein_ids[is.na(combined_hgt_phrog_raw$sister_clade_protein_ids)] <- ""
combined_hgt_phrog_raw$cousin_clade_protein_ids[is.na(combined_hgt_phrog_raw$cousin_clade_protein_ids)] <- ""
combined_hgt_phrog_raw$grandma_clade_protein_ids[is.na(combined_hgt_phrog_raw$grandma_clade_protein_ids)] <- ""

# Convert infinite values to 0
combined_hgt_phrog_raw$sister_clade_fraction_euk_species[is.infinite(combined_hgt_phrog_raw$sister_clade_fraction_euk_species)] <- 0
combined_hgt_phrog_raw$cousin_clade_fraction_euk_species[is.infinite(combined_hgt_phrog_raw$cousin_clade_fraction_euk_species)] <- 0

# Function to find largest parenthetical element in comma-separated string
get_largest_elems <- function(input_string, min_threshold) {
  if (is.na(input_string)) {
    return("")
  }
  if (input_string == "") {
    return("")
  }
  
  # Split the string into individual elements
  elements <- strsplit(input_string, ",")[[1]]
  
  # Remove empty elements
  elements <- elements[elements != ""]
  
  # Extract the names and the numeric values inside parentheses
  name_values <- gsub(".*\\(([^)]+)\\)", "\\1", elements)  # Extract numbers
  names <- gsub("\\(.*\\)", "", elements)  # Extract names
  
  # Convert the numeric values to a numeric vector
  num_values <- as.numeric(name_values)
  
  # Find the maximum number
  max_value <- max(num_values)
  
  # Extract the names corresponding to the maximum value
  max_names <- names[num_values == max_value]
  
  if (max_value >= min_threshold) {
    return(paste0(paste0(max_names, "(", max_value, ")"), collapse=","))
  } else {
    return("")
  }
}


## Identify HGT donors

# Set minimum relative fraction species coverage per prokaryotic group. Require enrichment equal to at least the baseline frequency
fraction_species_coverage_threshold <- 1

combined_hgt_phrog <- combined_hgt_phrog_raw

# Filter by presence of prokaryotes in sister or cousin as evidence of prok->euk HGT
combined_hgt_phrog <- combined_hgt_phrog %>% filter(sister_clade_protein_ids != "" | cousin_clade_protein_ids != "")

## Compute fraction of species per group, adjusted for number of species in group
# Only sister clade
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog %>% rowwise() %>% mutate(self_sister_clade_prok_species_in_clades = sister_clade_prok_species_in_clades)
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% separate_rows(self_sister_clade_prok_species_in_clades, sep=",") %>% filter(self_sister_clade_prok_species_in_clades != "") %>% group_by(PhROG_id) %>% filter(!duplicated(self_sister_clade_prok_species_in_clades))
combined_hgt_phrog_ancestor_prok_origin_species_counts$tree_id <- uniprot_proteomes_all_tax$tree_id[match(combined_hgt_phrog_ancestor_prok_origin_species_counts$self_sister_clade_prok_species_in_clades, uniprot_proteomes_all_tax$TaxId)]
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% group_by(PhROG_id, tree_id) %>% summarize(n_species = n())
uniprot_proteomes_all_tax_prokaryote_species_counts <- uniprot_proteomes_all_tax %>% filter(domain %in% c("Bacteria", "Archaea")) %>% group_by(tree_id) %>% summarize(total_species_count = n())
combined_hgt_phrog_ancestor_prok_origin_species_counts$total_species_count <- uniprot_proteomes_all_tax_prokaryote_species_counts$total_species_count[match(combined_hgt_phrog_ancestor_prok_origin_species_counts$tree_id, uniprot_proteomes_all_tax_prokaryote_species_counts$tree_id)]
# Fraction of clade corrected for prokgroup size
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% group_by(PhROG_id) %>% mutate(raw_fraction_species_coverage = (n_species / sum(n_species)), fraction_group_species_of_total_species = total_species_count / sum(uniprot_proteomes_all_tax_prokaryote_species_counts$total_species_count)) %>% mutate(fraction_species_coverage = raw_fraction_species_coverage / fraction_group_species_of_total_species)
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts[order(combined_hgt_phrog_ancestor_prok_origin_species_counts$fraction_species_coverage, decreasing=TRUE),]
combined_hgt_phrog_ancestor_prok_origin_species_counts <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% rowwise() %>% mutate(combined_label = paste0(tree_id, "(", fraction_species_coverage, ")"))
combined_hgt_phrog_ancestor_prok_origin_species_counts_summary <- combined_hgt_phrog_ancestor_prok_origin_species_counts %>% group_by(PhROG_id) %>% summarize(HGT_self_sister = paste0(combined_label, collapse=","))
combined_hgt_phrog$HGT_self_sister <- combined_hgt_phrog_ancestor_prok_origin_species_counts_summary$HGT_self_sister[match(combined_hgt_phrog$PhROG_id, combined_hgt_phrog_ancestor_prok_origin_species_counts_summary$PhROG_id)]
combined_hgt_phrog <- combined_hgt_phrog %>% filter(!is.na(HGT_self_sister)) %>% filter(HGT_self_sister != "")
# Assign most likely prokgroup by highest fraction species coverage
combined_hgt_phrog <- combined_hgt_phrog %>% rowwise() %>% mutate(HGT_self_sister = get_largest_elems(HGT_self_sister, fraction_species_coverage_threshold))

# Get largest parent PhROG@LECA per ACANB protein
combined_hgt_phrog_aca <- combined_hgt_phrog %>% rowwise() %>% filter(grepl("1257118_", reference_protein_ids, fixed=TRUE)) %>% rowwise() %>% mutate(protein_id = paste(reference_protein_ids, nonvertical_protein_ids, sep=",")) %>% separate_rows(protein_id, sep=",")

# Get just the ACANB proteins in their primary OGs and largest PhROGs (allow nonvertical proteins from other species that are secondary HGT)
selected_nodes <- c("1257118", "Node201_Amoebozoa", "Node141_Amorphea")
parent_progs_long_filter <- parent_progs_long %>% filter(label %in% selected_nodes)
parent_progs_long_filter <- parent_progs_long_filter %>% filter(grepl("1257118_", protein_id, fixed=TRUE))
combined_hgt_phrog_aca <- merge(combined_hgt_phrog_aca, parent_progs_long_filter, by.x=c("protein_id", "OG_id", "label"), by.y=c("protein_id", "OG_id", "label"))

# Filter by absence of eukaryotes in sister+cousin clade as evidence of recent HGT from prokaryotes
combined_hgt_phrog_aca <- combined_hgt_phrog_aca %>% filter(sister_clade_fraction_euk_species == 0 & cousin_clade_fraction_euk_species == 0)

# Remove duplicates (i.e. same ACANB protein mapped to multiple acquisition events) by taking highest prokaryotic enrichment, breaking ties by posterior clade support
combined_hgt_phrog_aca <- combined_hgt_phrog_aca %>% mutate(top_enrichment_odds_ratio = as.numeric(gsub("\\)", "", gsub(".*\\(", "", HGT_self_sister))))
combined_hgt_phrog_aca <- combined_hgt_phrog_aca %>% distinct() %>% group_by(protein_id) %>% filter(top_enrichment_odds_ratio == max(top_enrichment_odds_ratio)) %>% filter(count == max(count)) %>% ungroup() %>% filter(!duplicated(protein_id))

# Get mito localizations
gold_gene_accession_OG_id_df <- read.table(here("data/mito_orthogroups", "mito_proteins_experimental.and.mtDNA_2026.04.05.tsv"), sep="\t", header=TRUE)
combined_hgt_phrog_aca$in_mitocarta <- 0
combined_hgt_phrog_aca$in_mitocarta[combined_hgt_phrog_aca$protein_id %in% gold_gene_accession_OG_id_df$gene_accession] <- 1

combined_hgt_phrog_aca_summary <- combined_hgt_phrog_aca %>% group_by(label) %>% summarize(count = n())
combined_hgt_phrog_mito_summary <- combined_hgt_phrog_aca %>% filter(in_mitocarta == 1) %>% group_by(label) %>% summarize(count = n())

myPalette <- brewer.pal(3, "Set2")
pdf(here('Acanthamoeba_whole_proteome_HGT_pie.pdf'), width=8, height=4)
pie(combined_hgt_phrog_aca_summary$count, labels = paste0(combined_hgt_phrog_aca_summary$label, "(", combined_hgt_phrog_aca_summary$count, ")"), col=myPalette)
dev.off()

pdf(here('Acanthamoeba_mitoproteome_HGT_pie.pdf'), width=8, height=4)
pie(combined_hgt_phrog_mito_summary$count, labels = paste0(combined_hgt_phrog_mito_summary$label, "(", combined_hgt_phrog_mito_summary$count, ")"), col=myPalette)
dev.off()

selected_columns <- c("protein_id", "in_mitocarta", "label", "OG_id", "count", "reference_protein_ids", "HGT_self", "HGT_sister", "HGT_cousin", "HGT_self_sister")
combined_hgt_phrog_aca_out <- combined_hgt_phrog_aca[,selected_columns]
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "label"] <- "HGT_timing"
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "count"] <- "support"
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "reference_protein_ids"] <- "eukaryotic_homologs"
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "HGT_self"] <- "self_clade_prokaryote_groups(fraction_of_species_present)"
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "HGT_sister"] <- "sister_clade_prokaryote_groups(fraction_of_species_present)"
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "HGT_cousin"] <- "cousin_clade_prokaryote_groups(fraction_of_species_present)"
colnames(combined_hgt_phrog_aca_out)[colnames(combined_hgt_phrog_aca_out) == "HGT_self_sister"] <- "top_enriched_prokaryote_group(enrichment_odds_ratio)"
write.table(combined_hgt_phrog_aca_out, here('combined_HGT_results_Acanthamoeba.tsv'), sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

# Summarize HGT donors by prokgroup
combined_hgt_phrog_aca_donors <- combined_hgt_phrog_aca %>% separate_rows(HGT_self_sister, sep=",") %>% rowwise() %>% mutate(taxon = gsub("\\(.*\\)", "", HGT_self_sister), fraction_species_coverage = as.numeric(gsub(".*\\(([^)]+)\\)", "\\1", HGT_self_sister))) %>% filter(fraction_species_coverage >= fraction_species_coverage_threshold) %>% group_by(taxon, label) %>% summarize(count = length(unique(PhROG_id)), PhROG_id_collapse = paste0(PhROG_id, collapse=","))

combined_hgt_phrog_aca_donors_counts <- combined_hgt_phrog_aca_donors %>% group_by(taxon) %>% summarize(total_count = sum(count))
combined_hgt_phrog_aca_donors_counts <- combined_hgt_phrog_aca_donors_counts[order(combined_hgt_phrog_aca_donors_counts$total_count, decreasing=TRUE),]
combined_hgt_phrog_aca_donors$taxon <- factor(combined_hgt_phrog_aca_donors$taxon, levels = combined_hgt_phrog_aca_donors_counts$taxon)
combined_hgt_phrog_aca_donors$acquisition_timing <- combined_hgt_phrog_aca_donors$label
combined_hgt_phrog_aca_donors$acquisition_timing[combined_hgt_phrog_aca_donors$acquisition_timing == "1257118"] <- "Acanthamoeba"
combined_hgt_phrog_aca_donors$acquisition_timing[combined_hgt_phrog_aca_donors$acquisition_timing == "Node201_Amoebozoa"] <- "Amoebozoa"
combined_hgt_phrog_aca_donors$acquisition_timing[combined_hgt_phrog_aca_donors$acquisition_timing == "Node141_Amorphea"] <- "Amorphea"
combined_hgt_phrog_aca_donors$acquisition_timing <- factor(combined_hgt_phrog_aca_donors$acquisition_timing, levels = c("Acanthamoeba", "Amoebozoa", "Amorphea"))

pdf(here('Acanthamoeba_HGT_donors_bar.pdf'), width=8, height=6)
ggplot(data = combined_hgt_phrog_aca_donors, aes(x = taxon, y = count, fill = acquisition_timing)) + geom_col() + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1)) + ylab("Horizontal gene transfers") + coord_flip()
dev.off()

