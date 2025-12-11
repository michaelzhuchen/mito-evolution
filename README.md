# mito_evolution

## Overview

This repository contains code used in the study “Comparative mitochondrial proteomics across the tree of life reveals evolutionary, functional, and therapeutic insights”, part of the MitoCarta Tree of Life (MitoTOL) project, for orthogroup inference, phylogenetic profiling, ancestral reconstructions, eukaryogenesis timing, and comparative analyses of mitochondrial proteomes across the tree of life. This repository contains custom code and scripts used to generate the datasets released in the associated Zenodo dataset (```10.5281/zenodo.17823714```).

## Repository structure

```
mito-evolution/
├── abSENSE_HMM/ # estimate homology detection power across orthogroups
├── branch_length_timing/ # estimate origination timing from branch length
├── deeploc/ # retrain DeepLoc on new mitoproteomes and predict
├── horizontal_gene_transfer/ # identify prokaryote-derived HGT and putative donor lineages
├── orthogroups/ # refined orthogroups inference pipeline
├── phylogenetic_profiling/ # create OG/PhROG based homology matrix for phylogenetic profiling by CLIME
├── phylogenetically_resolved_orthogroups/ # phylogenetically-resolved orthogroup inference pipeline from reconciled trees
├── prokaryote_phenotype/ # prokaryote phenotype enrichment analysis
├── protein_phylogeny/ # protein phylogeny inference, processing, and reconciliation pipeline
├── species_phylogeny/ # species tree inference
├── ancestral_reconstruction/ # ancestral (mito)proteome reconstruction
├── utils/ # utility scripts
├── LICENSE # MIT License
└── README.md # this file
```


### Prerequisites

- Bash, Python (version ≥ 3.8), R (version ≥ 4.1)  
- Required packages and dependencies — see per-module documentation for details 
- Computational environment: To ensure efficiency and feasibility of computationally-intensive phylogenomic analysis, highly parallelizable scripts in the repository (ending in `_uger.sh`) are specifically designed to be run in parallel using taskarrays on an UGER-managed HPC, while other computationally-intensive and highly parallelizable scripts (ending in `_gcp.sh`) are specifically designed to be run as a Batch job on Google Cloud Compute Engine.


### Typical workflow
1. **Clone repository**
   ```bash
   git clone https://github.com/michaelzhuchen/mito-evolution.git
   ```

2. **Download data**  
   Download the full dataset from the Zenodo archive (```10.5281/zenodo.17823714```) and unpack each tar.gz archive to a directory with the same name within the `mito-evolution` directory, e.g.:

   ```bash
   tar -xvf data.tar.gz -C mito-evolution/data
   ```

3. **Install dependencies**  
   Install required packages and environments - see per-module documentation for details

4. **Run scripts**  
   Navigate to the module and run scripts, e.g.:
   ```bash
   cd mito-evolution/ancestral_reconstruction
   Rscript reconstruction.R
   ```



