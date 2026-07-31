# mito_evolution

## Overview

This repository contains code used in the study "Comparative analysis of mitochondrial proteomes across the tree of life”, part of the MitoCarta Tree of Life (MitoTOL) project, for orthogroup inference, ancestral reconstructions, eukaryogenesis timing, and comparative analyses of mitochondrial proteomes across the tree of life. This repository contains custom code and scripts used to generate the datasets released in the associated Zenodo dataset (```doi.org/10.5281/zenodo.21632671```).

## Repository structure

```
mito-evolution/
├── abSENSE_HMM/ # estimate homology detection power across orthogroups
├── benchmarks/ # benchmarks for orthogroups and phylogenetically-resolved orthogroups
├── branch_length_timing/ # estimate origination timing from branch length
├── comparative/ # comparative mitoproteome analyses
├── deeploc/ # train DeepLoc2.0-mito on new mitoproteomes and predict
├── horizontal_gene_transfer/ # identify prokaryote-derived HGT and putative donor lineages
├── orthogroups/ # refined orthogroups inference pipeline
├── phylogenetically_resolved_orthogroups/ # phylogenetically-resolved orthogroup inference pipeline from reconciled trees
├── prokaryote_phenotype/ # prokaryote phenotype enrichment analysis
├── protein_phylogeny/ # protein phylogeny inference, processing, and reconciliation pipeline
├── reconstruction/ # ancestral reconstruction
├── species_phylogeny/ # species tree inference
├── specific_analyses_per_organism/ # MitoTOL organism-specific analyses
├── utils/ # utility scripts
├── LICENSE # MIT License
└── README.md # this file
```


### Prerequisites

- R (version 4.3.1)
- Python (version 3.9.15)
- Required packages and dependencies — see per-module documentation for details


### Typical workflow
1. **Clone repository**
   ```bash
   git clone https://github.com/michaelzhuchen/mito-evolution.git
   ```
   Note that hereafter ```mito-evolution``` is used as shorthand for the file path of this base directory.

2. **Download data**  
   Download the full dataset from the Zenodo archive and unpack each tar.xz archive to a directory with the same name within the `mito-evolution` directory:

   ```bash
   tar -xJf alignments_and_initial_trees.tar.xz -C mito-evolution/alignments_and_initial_trees
   tar -xJf data.tar.xz -C mito-evolution/data
   tar -xJf pruned_rooted_trees.tar.xz -C mito-evolution/pruned_rooted_trees
   tar -xJf reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.1.tar.xz -C mito-evolution/reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.1
   tar -xJf reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.2.tar.xz -C mito-evolution/reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.2
   tar -xJf reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.3.tar.xz -C mito-evolution/reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.3
   tar -xJf reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.4.tar.xz -C mito-evolution/reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.4
   tar -xJf reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.5.tar.xz -C mito-evolution/reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.5
   tar -xJf reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.6.tar.xz -C mito-evolution/reconciled_consensus_trees_branch_length_optimization_with_supports_species.tree.6
   tar -xJf reconciled_consensus_trees_for_timing_species.tree.1.tar.xz -C mito-evolution/reconciled_consensus_trees_for_timing_species.tree.1
   tar -xJf reconciled_consensus_trees_for_timing_species.tree.2.tar.xz -C mito-evolution/reconciled_consensus_trees_for_timing_species.tree.2
   tar -xJf reconciled_consensus_trees_for_timing_species.tree.3.tar.xz -C mito-evolution/reconciled_consensus_trees_for_timing_species.tree.3
   tar -xJf reconciled_consensus_trees_for_timing_species.tree.4.tar.xz -C mito-evolution/reconciled_consensus_trees_for_timing_species.tree.4
   tar -xJf reconciled_consensus_trees_for_timing_species.tree.5.tar.xz -C mito-evolution/reconciled_consensus_trees_for_timing_species.tree.5
   tar -xJf reconciled_consensus_trees_for_timing_species.tree.6.tar.xz -C mito-evolution/reconciled_consensus_trees_for_timing_species.tree.6
   tar -xJf reconciled_trees_posterior_clades_species.tree.1.tar.xz -C mito-evolution/reconciled_trees_posterior_clades_species.tree.1
   tar -xJf reconciled_trees_posterior_clades_species.tree.2.tar.xz -C mito-evolution/reconciled_trees_posterior_clades_species.tree.2
   tar -xJf reconciled_trees_posterior_clades_species.tree.3.tar.xz -C mito-evolution/reconciled_trees_posterior_clades_species.tree.3
   tar -xJf reconciled_trees_posterior_clades_species.tree.4.tar.xz -C mito-evolution/reconciled_trees_posterior_clades_species.tree.4
   tar -xJf reconciled_trees_posterior_clades_species.tree.5.tar.xz -C mito-evolution/reconciled_trees_posterior_clades_species.tree.5
   tar -xJf reconciled_trees_posterior_clades_species.tree.6.tar.xz -C mito-evolution/reconciled_trees_posterior_clades_species.tree.6
   tar -xJf reconciled_trees_species.tree.1.tar.xz -C mito-evolution/reconciled_trees_species.tree.1
   tar -xJf reconciled_trees_species.tree.2.tar.xz -C mito-evolution/reconciled_trees_species.tree.2
   tar -xJf reconciled_trees_species.tree.3.tar.xz -C mito-evolution/reconciled_trees_species.tree.3
   tar -xJf reconciled_trees_species.tree.4.tar.xz -C mito-evolution/reconciled_trees_species.tree.4
   tar -xJf reconciled_trees_species.tree.5.tar.xz -C mito-evolution/reconciled_trees_species.tree.5
   tar -xJf reconciled_trees_species.tree.6.tar.xz -C mito-evolution/reconciled_trees_species.tree.6
   tar -xJf species_fastas.tar.xz -C mito-evolution/species_fastas
   ```

3. **Install dependencies**  
   Install required packages and environments - see per-module documentation for details.

4. **Run scripts**  
   Follow the per-module documentation where provided, or navigate to the module and run scripts from that working directory, e.g.:
   ```bash
   cd mito-evolution/ancestral_reconstruction
   Rscript reconstruction.R
   ```

   Caution: some scripts generate output files that may overwrite preexisting files downloaded from Zenodo.


