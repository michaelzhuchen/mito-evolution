# benchmarks

## Overview

This folder contains scripts used to perform benchmark analysis.


### Typical workflow
1. **Install dependencies**  
   Install required packages and environments.

   Install R packages within your R environment:
   ```R
   install.packages(here)
   install.packages(tidyverse)
   install.packages(ggrepel)
   install.packages(gplots)
   install.packages(basetheme)
   install.packages(RColorBrewer)
   install.packages(reshape2)
   install.packages(ape)
   install.packages(castor)
   install.packages(phytools)
   ```


2. **Run scripts**  
   Navigate to the module and run scripts from the working directory, e.g.:
   ```bash
   cd mito-evolution/benchmarks

   # Run scripts
   Rscript orthogroups_cryoEM_pathways_benchmark.R # Run cryoEM complexes/pathways benchmark

   Rscript PhROG_HGT_parameter_control.R # Run HGT inference for PhROGs benchmark
   
   # Run duplication inference for PhROGs benchmark over a range of species overlap thresholds
   Rscript PhROG_duplication_parameter_control.R species.overlap.0.0
   Rscript PhROG_duplication_parameter_control.R species.overlap.0.01
   Rscript PhROG_duplication_parameter_control.R species.overlap.0.1
   Rscript PhROG_duplication_parameter_control.R species.overlap.0.2
   Rscript PhROG_duplication_parameter_control.R species.overlap.0.3

   ```


