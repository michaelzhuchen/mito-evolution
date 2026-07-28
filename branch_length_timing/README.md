# branch_length_timing

## Overview

This folder contains scripts used to perform branch length timing analysis.


### Typical workflow
1. **Install dependencies**  
   Install required packages and environments.

   Install R packages within your R environment:
   ```R
   install.packages(here)
   install.packages(tidyverse)
   install.packages(gplots)
   install.packages(basetheme)
   install.packages(RColorBrewer)
   install.packages(ggtree)
   install.packages(ape)
   install.packages(castor)
   install.packages(phytools)
   ```


2. **Run scripts**  
   Navigate to the module and run scripts from the working directory, e.g.:
   ```bash
   cd mito-evolution/branch_length_timing

   # Run scripts
   DATASET_NAME="species_tree_1" # Select species tree dataset, e.g. species_tree_[1-6]
   SELECTED_TAX_LEVEL="Node34_Eukaryota" # Select focal species tree node, e.g. Node34_Eukaryota, Node39_Archaeplastida
   OGID="MOG0001047" # Select orthogroup

   Rscript extract_normalized_branch_lengths_parallel.R $DATASET_NAME $SELECTED_TAX_LEVEL $OGID # Extract branch lengths from reconciled consensus trees with optimized branch lengths

   Rscript timing_from_branch_lengths.R $DATASET_NAME # Analyze branch length distributions across orthogroups

   ```


