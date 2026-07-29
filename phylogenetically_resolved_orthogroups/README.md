# phylogenetically_resolved_orthogroups

## Overview

This folder contains scripts used to perform phylogenetically resolved orthogroup (PhROG) analysis.


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
   cd mito-evolution/phylogenetically_resolved_orthogroups

   # Get posterior distribution of clades in protein phylogenies
   DATASETNAME="species_tree_1"
   RECTREESDIR="reconciled_trees_species.tree.1"
   RECTREESPOSTERIORDIR="reconciled_trees_posterior_clades_species.tree.1"
   OGID="MOG0001047"
   Rscript get_posterior_clades.R $DATASETNAME $RECTREESDIR $RECTREESPOSTERIORDIR $OGID

   # Infer PhROGs
   Rscript generate_PhROGs.R $DATASETNAME $RECTREESPOSTERIORDIR $OGID

   # Attach ancestral mitochondrial localization to PhROGs
   RECTREESCONSENSUSDIR="reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.1"
   Rscript infer_mito_localization_for_PhROGs.R $DATASETNAME $RECTREESDIR $RECTREESCONSENSUSDIR $OGID

   ```


