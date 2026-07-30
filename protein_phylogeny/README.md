# protein_phylogeny

## Overview

This folder contains scripts used to perform protein phylogenetic inference.


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

   Install muscle (v5), available here: https://www.drive5.com/muscle/. \
   Install IQ-TREE (v2.3.6), available here: https://iqtree.github.io/. \

2. **Run scripts**  
   Navigate to the base directory and run scripts. Note that these steps can be compute-intensive and may run slowly on large datasets.
   ```bash
   cd mito-evolution

   OGID="MOG0001047" # Select orthogroup
   DATASETNAME="species_tree_1"
   N_CORES="2"

   # Multiple sequence alignment
   chmod u+x protein_phylogeny/muscle_trim.sh
   ./protein_phylogeny/muscle_trim.sh $OGID $N_CORES

   # Maximum-likelihood tree inference with sample of BNNI-optimized UFBOOT trees
   iqtree2 -s ${OGID}.faa_clipkit.gappy.msa -m MFP -cmin 4 -cmax 12 -nt $N_CORES -bb 1000 -alrt 1000 -seed 42 -wbtl -bnni

   # Infer origin, root, and remove potential prokaryotic contaminants
   Rscript protein_phylogeny/infer_gene_family_origin_prune_root_trees.R ${OGID}.faa_clipkit.gappy.msa $DATASETNAME

   # Reconcile with species tree
   ALERAXDIR="reconciled_trees_species.tree.1"
   chmod u+x protein_phylogeny/run_alerax.sh
   ./protein_phylogeny/run_alerax.sh $OGID $DATASETNAME $ALERAXDIR

   # Infer branch lengths for majority-rule consensus tree
   CONSENSUSBLOPTDIR="reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.1"
   chmod u+x protein_phylogeny/optimize_consensus_branchlength.sh
   ./protein_phylogeny/optimize_consensus_branchlength.sh $OGID $ALERAXDIR $CONSENSUSBLOPTDIR
   ```


