# species_phylogeny

## Overview

This folder contains scripts used to infer the species phylogeny.


### Typical workflow
1. **Install dependencies**  
   Install required packages and environments.

   Install R packages within your R environment:
   ```R
   install.packages(here)
   install.packages(tidyverse)
   install.packages(Biostrings)
   install.packages(basetheme)
   install.packages(RColorBrewer)
   install.packages(ape)
   install.packages(phytools)
   install.packages(TreeTools)
   ```

   Install IQ-TREE (v2.3.6) from https://iqtree.github.io/.
   

2. **Run scripts**  
   Navigate to the module and run scripts from the working directory, e.g.:
   ```bash
   cd mito-evolution/species_phylogeny

   # Prepare concatenated MSA
   Rscript concatenate_msa.R

   # Infer species tree. Note that this can take some time, so multiple cores are recommended.
   chmod u+x infer_species_tree_iqtree2.sh
   ./infer_species_tree_iqtree2.sh

   # Process species tree
   Rscript process_species_tree.R

   ```


