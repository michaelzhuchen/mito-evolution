# mito_evolution

## Overview

This repository contains code used in the study "Comparative analysis of mitochondrial proteomes across the tree of life”, part of the MitoCarta Tree of Life (MitoTOL) project, for orthogroup inference, ancestral reconstructions, eukaryogenesis timing, and comparative analyses of mitochondrial proteomes across the tree of life. This repository contains custom code and scripts used to generate the datasets released in the associated Zenodo dataset (```doi.org/10.5281/zenodo.20672825```).

## Repository structure

```
mito-evolution/
├── abSENSE_HMM/ # estimate homology detection power across orthogroups
├── branch_length_timing/ # estimate origination timing from branch length
├── deeploc/ # train DeepLoc2.0-mito on new mitoproteomes and predict
├── horizontal_gene_transfer/ # identify prokaryote-derived HGT and putative donor lineages
├── orthogroups/ # refined orthogroups inference pipeline
├── phylogenetically_resolved_orthogroups/ # phylogenetically-resolved orthogroup inference pipeline from reconciled trees
├── prokaryote_phenotype/ # prokaryote phenotype enrichment analysis
├── protein_phylogeny/ # protein phylogeny inference, processing, and reconciliation pipeline
├── reconstruction/ # ancestral reconstruction
├── species_phylogeny/ # species tree inference
├── utils/ # utility scripts
├── LICENSE # MIT License
└── README.md # this file
```


### Prerequisites

- Bash, Python (version ≥ 3.8), R (version ≥ 4.1)  
- Required packages and dependencies — see per-module documentation for details


### Typical workflow
1. **Clone repository**
   ```bash
   git clone https://github.com/michaelzhuchen/mito-evolution.git
   ```

2. **Download data**  
   Download the full dataset from the Zenodo archive (```doi.org/10.5281/zenodo.20672825```) and unpack each tar.xz archive to a directory with the same name within the `mito-evolution` directory, e.g.:

   ```bash
   tar -xJf data.tar.xz -C mito-evolution/data
   ```

3. **Install dependencies**  
   Install required packages and environments - see per-module documentation for details

4. **Run scripts**  
   Navigate to the module and run scripts from that working directory, e.g.:
   ```bash
   cd mito-evolution/ancestral_reconstruction
   Rscript reconstruction.R
   ```



