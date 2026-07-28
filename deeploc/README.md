# deeploc

## Overview

This folder contains scripts used to train and generate predictions using DeepLoc2.0-mito.


### Typical prediction workflow
1. **Install dependencies**  

   Download DeepLoc Version 2.0 tar archive from https://services.healthtech.dtu.dk/services/DeepLoc-2.0/ to a desired installation location.

   Extract the downloaded tar archive. Note that hereafter ```deeploc2_package``` is used as shorthand for the full file path of the extracted directory:
   ```bash
   tar -xvf deeploc-2.0.All.tar.gz
   ```
   
   Install virtual environment (if not already installed with python installation).

   Use python v3.9.15 to create a virtual environment .venv within deeploc2_package:
   ```bash
   cd deeploc2_package
   python -m venv .venv
   ```

   Activate virtual environment and install dependencies: 
   ```bash
   source .venv/bin/activate

   # Install numpy v1.26.4
   python -m pip install numpy==1.26.4

   # Install sentencepiece v0.2.0
   python -m pip install sentencepiece==0.2.0

   # Install legacy setuptools version<82 with pkg_resources 
   python -m pip install "setuptools<82"
   ```

   Extract data.tar.xz archive downloaded from Zenodo repository.
   ```bash
   tar -xJf data.tar.xz -C mito-evolution/data
   ```

2) Replace the default DeepLoc2.0 model weights with the weights from the fully-retrained DeepLoc2.0-mito model:
   ```bash
   rm -r DeepLoc2/models/models_prott5
   cp -r mito-evolution/data/deeploc/models/models_prott5_PCP_All_2026.04.05_graphpart.expect1.nopriority DeepLoc2/models/models_prott5
   ```

3) Compile DeepLoc2.0 using the updated model weights:
   ```bash
   rm -rf build
   python -m pip install .
   ```
   
2. **Run scripts**  
   Generate predictions using the retrained DeepLoc2.0 model. Note that the first run will take longer since protein language model files will need to be downloaded (cached thereafter).
   ```bash
   cd deeploc2_package
   source .venv/bin/activate
   deeploc2 -f test.fasta -o outputs -m Accurate
   ```


### Typical training workflow
1. **Install dependencies**  
   Install required packages and environments.

   Install FASTA package (v36.3.8i). Precompiled binaries available at https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36/, or install from source at https://github.com/wrpearson/fasta36/releases/tag/v36.3.8i_14-Nov-2020.

   Install R packages within your R environment:
   ```R
   install.packages(here)
   install.packages(tidyverse)
   install.packages(Biostrings)
   ```

   Extract data.tar.xz archive downloaded from Zenodo repository.
   ```bash
   tar -xJf data.tar.xz -C mito-evolution/data
   ```

   Install Anaconda. We recommend the Miniconda distribution, available here: https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install

   Install DeepLoc2.0 development version from https://github.com/teevee112/DeepLoc-2.0:
   ``` bash
   git clone https://github.com/teevee112/DeepLoc-2.0.git
   cd DeepLoc-2.0
   conda env create -f environment.yml
   ```

   Replace original scripts:
   ```bash
      cp mito-evolution/deeploc/2_train/data.py src/data.py
      cp mito-evolution/deeploc/2_train/metrics.py src/metrics.py
      cp mito-evolution/deeploc/2_train/model.py src/model.py
   ```
   
   Modify updated scripts:
      ```src/constants.py```: update filepaths to point to new training dataset
      ```data_files/embed_configs/swissprot_t5.yaml```: update filepaths to point to new training dataset
      ```src/model.py```: update ```pos_weights_bce``` vector (for localization label weights) in ```focal_loss()``` function based on the frequency of each label in the training dataset

2. **Run scripts**  
   
   Train DeepLoc2.0:
   ```bash
   conda activate deeploc20

   # Train subcellular localization module
   python train_sl.py --model Accurate

   # Train sorting signal module
   python train_ss.py --model Accurate
   ```

   After training, the trained model weights are available in ```models/```.




