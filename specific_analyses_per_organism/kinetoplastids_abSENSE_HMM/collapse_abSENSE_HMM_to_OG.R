library(tidyverse)

dataset_name <- "species_tree_1"

# Collapse abSENSE-HMM results to single label per OG for core species
progs_long <- read.table(here("data/phylogenetically_resolved_orthogroups", dataset_name, "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(progs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")

progs_long <- progs_long %>% mutate(Orthogroup = PROG_id, OG_id = gsub("_Node.*", "", PROG_id), taxid = gsub("_.*", "", protein_id))
PhROG_OGs <- unique(progs_long$OG_id)

# Read in OGs
ogs_long <- read.table(here("data/orthogroups/refined_orthogroups", "refined_OGs_euk203spp_long.txt"), sep="\t", header=TRUE)
colnames(ogs_long) <- c("accession", "Orthogroup", "taxid", "BOOL_PRIMARY_OG")
core_species_taxids <- c("9606", "559292", "1257118", "5741", "185431", "5689", "3702", "32595")
ogs_long_core_species <- ogs_long %>% filter(taxid %in% core_species_taxids)
core_species_OG_ids <- unique(ogs_long_core_species$Orthogroup)

# Results from mito OGs and OGs selected for tree-building
homology_power_agg_PhROG <- read.table(here("data/abSENSE_HMM", dataset_name, "absense_results.tsv"), sep="\t", header=TRUE)
colnames(homology_power_agg_PhROG) <- c("PhROG_id", "human_accessions", "L", "R", "OG_mrca_label", "OG_outgroup_mrca_label", "outgroup_kingdom_species", "fraction_of_detectable_outgroup_kingdoms", "probability_of_detection_in_any_outgroup_species", "hmm_length", "r_squared", "comments")
homology_power_agg_PhROG <- homology_power_agg_PhROG %>% rowwise() %>% mutate(OG_id = gsub("_Node.*", "", PhROG_id))

# Results from OGs not selected for tree-building
homology_power_agg_OG <- read.table(here("specific_analyses_per_organism", "kinetoplastids_abSENSE_HMM", "absense_OG_core_species_results.tsv"), sep="\t", header=FALSE)
colnames(homology_power_agg_OG) <- c("PhROG_id", "human_accessions", "L", "R", "OG_mrca_label", "OG_outgroup_mrca_label", "outgroup_kingdom_species", "fraction_of_detectable_outgroup_kingdoms", "probability_of_detection_in_any_outgroup_species", "hmm_length", "r_squared", "comments")
homology_power_agg_OG <- homology_power_agg_OG %>% rowwise() %>% mutate(OG_id = gsub("_Node.*", "", PhROG_id))

# Remove duplicates and combine
homology_power_agg_OG <- homology_power_agg_OG %>% filter(!OG_id %in% PhROG_OGs)
homology_power_agg <- rbind(homology_power_agg_PhROG, homology_power_agg_OG)

# Power to detect in any outgroup kingdoms per OG
all_OG_ids_detectable <- unique(homology_power_agg$OG_id[which(homology_power_agg$fraction_of_detectable_outgroup_kingdoms > 0 | homology_power_agg$OG_mrca_label == "Node34_Eukaryota")])
OG_ids_detectable <- core_species_OG_ids[core_species_OG_ids %in% all_OG_ids_detectable]
OG_ids_undetectable <- core_species_OG_ids[!core_species_OG_ids %in% all_OG_ids_detectable]

## Write out absense table, collapsed per OG
absense_OG_table <- rbind(data.frame(OG_id = OG_ids_detectable, label = "Powered"), data.frame(OG_id = OG_ids_undetectable, label = "Not powered"))
absense_OG_table <- absense_OG_table[order(absense_OG_table$OG_id),]
# write.table(absense_OG_table, here("absense_OG_core_species_results_collapsetoOG.tsv"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

