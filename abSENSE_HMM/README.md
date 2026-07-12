# mito_evolution

## Overview

This folder contains scripts used to perform abSENSE-HMM analysis. The analysis takes as input the phylogenetically-resolved orthogroups from the Zenodo archive (```doi.org/10.5281/zenodo.20672825```).


### Typical workflow
1. **Install dependencies**  
   Install required packages and environments.
   ```R (v4.3.1)```: https://cran.r-project.org/
   ```HMMER (v3.3.2)```: http://hmmer.org/documentation.html. Ensure that HMMER executables are in your path.

   Install R packages within the R environment:
   ```R
   install.packages(tidyverse)
   install.packages(Biostrings)
   install.packages(ape)
   install.packages(castor)
   ```


2. **Run scripts**  
   Navigate to the module and run scripts from the working directory, e.g.:
   ```bash
   cd mito-evolution/abSENSE_HMM
   
   # Set permissions
   chmod u+x hmm_LOO_PhROG.sh
   chmod u+x homology_power_PhROG.sh

   # Run scripts
   OGID="MOG0001047" # Select orthogroup
   ./hmm_LOO_PhROG.sh $OGID # Run leave-one-out hmmsearch
   ./homology_power_PhROG.sh $OGID # Run abSENSE-HMM model
   ```



