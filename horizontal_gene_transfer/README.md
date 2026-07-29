# horizontal_gene_transfer

## Overview

This folder contains scripts used to perform horizontal gene transfer (HGT) analysis.


### Typical workflow
1. **Install dependencies**  
   Install required packages and environments.

   Install R packages within your R environment:
   ```R
   install.packages(here)
   install.packages(tidyverse)
   install.packages(ape)
   install.packages(castor)
   install.packages(phytools)
   ```


2. **Run scripts**  
   Navigate to the module and run scripts from the working directory, e.g.:
   ```bash
   cd mito-evolution/horizontal_gene_transfer

   # Infer HGTs from posterior distribution of clades on reconciled phylogenies
   DATASETNAME="species_tree_1"
   INPUTDATA="reconciled_trees_posterior_clades_species.tree.1"
   SELECTEDTAXLEVEL="Node34_Eukaryota_parent" # for HGT gains into Eukaryota
   OGID="MOG0001047"
   Rscript get_HGT_from_posterior_clades.R $DATASETNAME $INPUTDATA $SELECTEDTAXLEVEL $OGID

   # Summarize HGT inferences
   DATASETNAME="species_tree_1"
   SELECTEDTAXLEVEL="Node34_Eukaryota" # include HGT gains at/above this ancestral species tree node.
   Rscript HGT_donors.R $DATASETNAME $SELECTEDTAXLEVEL

   ```
