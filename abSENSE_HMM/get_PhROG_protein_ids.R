suppressMessages(library(tidyverse))

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
phrogs_dir <- args[1]
out_dir <- args[2]
OG_id <- args[3]

phrogs_long <- read.table(here("data/phylogenetically_resolved_orthogroups", "PhROGs_long", "PhROGs_at_Node34_Eukaryota_parent_long.tsv"), sep="\t", header=TRUE)
colnames(phrogs_long) <- c("OG_id", "PROG_id", "label", "mito_localization_prob_mk", "mito_localization_prob_parsimony", "protein_id",  "BOOL_NONVERTICAL", "BOOL_primary_OG")

# Remove duplicates to assign proteins to their largest PhROG (ignore non-vertical proteins that are in a larger PhROG)
phrogs_long <- phrogs_long %>% filter(!duplicated(protein_id)) %>% filter(!grepl("dupelabel", protein_id, fixed=TRUE))
# Remove nonvertical proteins to only retain orthologs + in-paralogs
phrogs_long <- phrogs_long %>% filter(!BOOL_NONVERTICAL)
# Keep PhROGs that have at least 3 proteins (required for model fitting)
phrogs_long <- phrogs_long %>% group_by(PROG_id) %>% mutate(n_proteins = n()) %>% filter(n_proteins >= 3)

for (PROG_id_curr in unique(phrogs_long$PROG_id)) {
  phrogs_long_curr <- phrogs_long %>% filter(PROG_id == PROG_id_curr)
  write.table(phrogs_long_curr$protein_id, paste0(out_dir, "/", PROG_id_curr, "_LOO_protein_id_list.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
}

