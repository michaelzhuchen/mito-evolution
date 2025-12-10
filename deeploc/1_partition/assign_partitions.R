### Assign partitions based on graphpart output

# Load libraries
library(tidyverse)
library(Biostrings)


## Select training dataset
# all species
suffix <- "mitoPCP.mitoonlyPCP.nonmitoPCP_all"

# leave out Aca
# suffix <- "mitoPCP.mitoonlyPCP.nonmitoPCP_LOO.Aca"

# leave out kinetoplastids
# suffix <- "mitoPCP.mitoonlyPCP.nonmitoPCP_LOO.lta.tbr"

# leave out Ath
# suffix <- "mitoPCP.mitoonlyPCP.nonmitoPCP_LOO.ath"

# leave out Bdi
# suffix <- "mitoPCP.mitoonlyPCP.nonmitoPCP_LOO.bdi"

# leave out Giardia
# suffix <- "mitoPCP.mitoonlyPCP.nonmitoPCP_LOO.gla"


## Read in data
# Read in training dataset
multisub_labels_swissprot_mitoepi <- read.csv(here("data/deeploc/data_files", paste0("multisub_swissprot.retained_added.", suffix, ".csv")))

# Read in partitions from graphpart
graphpart_assignments <- read.csv(here("data/deeploc/graphpart", paste0("swissprot_", suffix, ".fasta_pident0.3_expect1_nopriority/graphpart_assignments.csv")))

# Assign partition
multisub_labels_swissprot_mitoepi$Partition <- graphpart_assignments$cluster[match(multisub_labels_swissprot_mitoepi$ACC, graphpart_assignments$AC)]
table(multisub_labels_swissprot_mitoepi$Partition)

# Randomize order
multisub_labels_swissprot_mitoepi$X <- sample(1:nrow(multisub_labels_swissprot_mitoepi))
multisub_labels_swissprot_mitoepi <- multisub_labels_swissprot_mitoepi[order(multisub_labels_swissprot_mitoepi$X),]

# Remove proteins without a partition (excluded due to cross-partition homology)
multisub_labels_swissprot_mitoepi_nomissing <- multisub_labels_swissprot_mitoepi %>% filter(!is.na(Partition))

# Set column names
colnames(multisub_labels_swissprot_mitoepi_nomissing) <- c("X", "ACC", "Kingdom", "Partition", "Membrane","Cytoplasm","Nucleus","Extracellular","Cell membrane","Mitochondrion","Plastid","Endoplasmic reticulum","Lysosome/Vacuole","Golgi apparatus","Peroxisome","Sequence","species","taxid")

# Fill missing values with nan
multisub_labels_swissprot_mitoepi_nomissing_nan <- multisub_labels_swissprot_mitoepi_nomissing
multisub_labels_swissprot_mitoepi_nomissing_nan[is.na(multisub_labels_swissprot_mitoepi_nomissing_nan)] <- "nan"

## Write out
# write.csv(multisub_labels_swissprot_mitoepi_nomissing_nan, here("data/deeploc/data_files", paste0("multisub_5partitions_swissprot.retained_added.", suffix, "_graphpart.expect1.nopriority_nomissing_nan.csv")), row.names=FALSE, quote=FALSE)


## Retrieve fasta for partitioned training dataset
# Read in sequences
swissprot_fasta <- readAAStringSet(here("data/deeploc/swissprot", "deeploc_swissprot.fasta"))
mitoepi_fasta <- readAAStringSet(here("data/deeploc/mitotol", "mitoepi_species_combined.fasta")) # all MITO-EPI species proteins
combined_fasta <- c(mitoepi_fasta, swissprot_fasta)
names(combined_fasta)[grep("^1257118_", names(combined_fasta))] <- gsub("_t.*", "", names(combined_fasta)[grep("^1257118_", names(combined_fasta))])

combined_fasta_for_train <- combined_fasta[multisub_labels_swissprot_mitoepi_nomissing$ACC]

# Trim sequences to maximum input length for ProtT5 model (4000 residues)
library(Biostrings)
maxLen  <- 4000L
halfLen <- maxLen %/% 2L

trim_middle <- function(x) {
  n <- nchar(as.character(x))
  if (n <= maxLen) {
    # leave short sequences untouched
    return(x)
  }
  # take first halfLen and last halfLen
  left  <- subseq(x, start=1L,                width=halfLen)
  right <- subseq(x, start=n - halfLen + 1L,  width=halfLen)
  xscat(left, right)
}

#    lapply gives a list of AAString objects, which we then re-wrap into an AAStringSet
trimmed_list <- lapply(combined_fasta_for_train, trim_middle)
combined_fasta_for_train_trimmed      <- AAStringSet(trimmed_list)
names(combined_fasta_for_train_trimmed) <- names(combined_fasta_for_train)

# writeXStringSet(combined_fasta_for_train_trimmed, here("data/deeploc/data_files", paste0("swissprot_", suffix, "_graphpart.expect1.nopriority_clipped4k.fasta")))


