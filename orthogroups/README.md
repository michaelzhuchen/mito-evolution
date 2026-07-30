# orthogroups

## Overview

This folder contains scripts used to perform orthogroup analysis. Note that this pipeline is compute-intensive for large numbers of species/proteins; many steps require parallel compute to finish in a reasonable timeframe.


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

   Install OrthoFinder (v2.5.5), using conda (recommended) or from source at https://github.com/davidemms/OrthoFinder.
   ```bash
   conda create --name OrthoFinder_v255
   conda activate OrthoFinder_v255
   conda install bioconda::orthofinder==2.5.5
   ```

   Install mmseqs2, available here: https://github.com/soedinglab/mmseqs2
   Install mafft, available here: https://mafft.cbrc.jp/alignment/software/
   Install ClipKIT, available here: https://github.com/jlsteenwyk/clipkit
   Install HH-suite, available here: https://github.com/soedinglab/hh-suite
   Install Open MPI, available here: https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html


   Extract tar archives downloaded from Zenodo repository.
   ```bash
   tar -xJf species_fastas.tar.xz -C mito-evolution/species_fastas
   tar -xJf data.tar.xz -C mito-evolution/data
   ```

2. **Run scripts**  
   Navigate to the base directory and run scripts, e.g.:
   ```bash
   cd mito-evolution
   ```

   # Raw orthogroups from OrthoFinder
   Run OrthoFinder to infer raw orthogroups from a set of species proteome FASTA files. Note that this is compute and memory intensive, so multiple CPUs and high RAM are recommended. The most time-intensive step is pairwise alignment, which can be precomputed and parallelized using the -op and -b options (see documentation here for details: https://github.com/davidemms/OrthoFinder)
   ```bash
   conda activate OrthoFinder_v255
   orthofinder -f species_fastas/ -t 256 -a 32 -og
   ```

   # Initially merged orthogroups from first merge
   Build alignments for each raw nonsingleton orthogroup, using FASTA files from OrthoFinder raw orthogroups as input. For each orthogroup, run:
   ```bash
   BASEDIR="first_merge"
   mkdir -p $BASEDIR
   RAWOGDIR='$BASEDIR/OG_raw'
   TMPDIR="tmp"
   OGID="OG0002514"

   # Downsample orthogroup using MMseqs2
   mmseqs easy-cluster $RAWOGDIR/${OGID}.faa $RAWOGDIR/${OGID} $TMPDIR --threads 1 --split-memory-limit 2G --min-seq-id 0.5 --cov-mode 2 -c 0.8 -v 1

   # Align using mafft
   mafft --auto --quiet --thread 1 --anysymbol "$RAWOGDIR/${OGID}_rep_seq.fasta" > "$RAWOGDIR/${OGID}.faa.msa"

   # Trim alignment using Clipkit
   conda activate clipkit_2.3.0
   clipkit -q -m gappy "$RAWOGDIR/${OGID}.faa.msa" -o "$RAWOGDIR/${OGID}.faa_clipkit.gappy.msa"
   python utils/remove_gap_only_seqs.py "$RAWOGDIR/${OGID}.faa_clipkit.gappy.msa"
   ```

   Build HH-suite HMM database from all raw nonsingleton orthogroup alignments.
   ```bash
   # Create folder and copy alignments as input for database construction
   DBNAME="OG_raw_trimmed_hhsuite"
   mkdir -p $BASEDIR/$DBNAME
   cp $RAWOGDIR/*_clipkit.gappy.msa $BASEDIR/$DBNAME

   # Build HH-suite database. Recommend multithreading with 8 cores using MPI.
   chmod u+x orthogroups/1_initially_merged_orthogroups/hhsearch/build_database_hhsuite.sh
   ./orthogroups/1_initially_merged_orthogroups/hhsearch/build_database_hhsuite.sh
   ```

   Run all-vs-all HMM alignments. For each orthogroup, run:
   ```bash
   OGID="OG0002514"
   QUERYFILE_SUFFIX=".faa_clipkit.gappy.msa"
   QUERYDIR="$BASEDIR/$DBNAME"
   TARGETDB="$BASEDIR/${DBNAME}_db"
   OUTDIR="$BASEDIR/hhsearch_out"

   chmod u+x 1_initially_merged_orthogroups/hhsearch/run_hhsearch.sh
   ./orthogroups/1_initially_merged_orthogroups/hhsearch/run_hhsearch.sh $OGID $QUERYFILE_SUFFIX $QUERYDIR $TARGETDB $OUTDIR

   # Compute overlap coefficient for filtered hits
   SUFFIX=".faa_clipkit.gappy.msa_global.hhr_expect1e-3"
   BOOL_STRUCTURE_MERGE="FALSE"
   Rscript orthogroups/helpers/compute_overlap_coef_for_tsv.R "$OUTDIR" "$OUTDIR" "$OGID" "$SUFFIX" "$BOOL_STRUCTURE_MERGE"
   ```

   Perform first merge and output the resulting initially merged orthogroups.
   ```bash
   Rscript orthogroups/1_initially_merged_orthogroups/merge/perform_first_merge.R
   Rscript orthogroups/1_initially_merged_orthogroups/merge/process_first_merged_orthogroups.R
   Rscript orthogroups/1_initially_merged_orthogroups/merge/get_initially_merged_orthogroups.R
   ```

   # Fully merged orthogroups from second merge
   
   Build alignments for each initially merged orthogroup, using FASTA files for initially merged orthogroups as input. For each orthogroup, run:
   ```bash
   BASEDIR="second_merge"
   mkdir -p $BASEDIR
   MERGEDOGDIR='$BASEDIR/OG_merged'
   TMPDIR="tmp"
   OGID="MOG0030340"

   # Align using mafft
   mafft --auto --quiet --thread 1 --anysymbol "$MERGEDOGDIR/${OGID}.faa" > "$MERGEDOGDIR/${OGID}.faa.msa"

   # Trim alignment using Clipkit
   conda activate clipkit_2.3.0
   clipkit -q -m gappy "$MERGEDOGDIR/${OGID}.faa.msa" -o "$MERGEDOGDIR/${OGID}.faa_clipkit.gappy.msa"
   python utils/remove_gap_only_seqs.py "$MERGEDOGDIR/${OGID}.faa_clipkit.gappy.msa"
   ```

   Build HH-suite HMM database from all initially merged orthogroup alignments.
   ```bash
   # Create folder and copy alignments as input for database construction
   DBNAME="OG_merged_trimmed_hhsuite"
   mkdir -p $BASEDIR/$DBNAME
   cp $RAWOGDIR/*_clipkit.gappy.msa $BASEDIR/$DBNAME

   # Build HH-suite database. Recommend multithreading with 8 cores using MPI.
   chmod u+x orthogroups/2_fully_merged_orthogroups/hhsearch/build_database_hhsuite.sh
   ./orthogroups/2_fully_merged_orthogroups/hhsearch/build_database_hhsuite.sh
   ```

   Run all-vs-all HMM alignments. For each orthogroup, run:
   ```bash
   OGID="MOG0030340"
   QUERYFILE_SUFFIX=".faa_clipkit.gappy.msa"
   QUERYDIR="$BASEDIR/$DBNAME"
   TARGETDB="$BASEDIR/${DBNAME}_db"
   OUTDIR="$BASEDIR/hhsearch_out"

   chmod u+x orthogroups/2_fully_merged_orthogroups/hhsearch/run_hhsearch.sh
   ./orthogroups/2_fully_merged_orthogroups/hhsearch/run_hhsearch.sh $OGID $QUERYFILE_SUFFIX $QUERYDIR $TARGETDB $OUTDIR

   # Compute overlap coefficient for filtered hits
   SUFFIX=".faa_clipkit.gappy.msa_global.hhr_prob30"
   BOOL_STRUCTURE_MERGE="TRUE"
   Rscript orthogroups/helpers/compute_overlap_coef_for_tsv.R "$OUTDIR" "$OUTDIR" "$OGID" "$SUFFIX" "$BOOL_STRUCTURE_MERGE"
   ```

   Run Foldseek alignments for each pair of orthogroups in the set of filtered HMM hits. First, obtain the predicted structures (e.g. from AFDB, Uniprot, or de novo prediction) for each pair of orthogroups in the set of filtered HMM hits to folders ```PDBDIR/[QUERY_OG]/structs_raw``` and ```PDBDIR/[TARGET_OG]/structs_raw```, corresponding to query and target orthogroups. Then, run the following scripts for each pair of orthogroups in the set of filtered HMM hits:
   ```bash
   QUERY_OG="MOG0030340"
   TARGET_OG="OG0059097"
   PDBDIR="PDBDIR"
   QUERYDIR="${QUERY_OG}/structs_raw"
   TARGETDIR="${TARGET_OG}/structs_raw"
   OUTDIR="$BASEDIR/foldseek_out"
   FSOUTFILE="$OUTDIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10.tsv"
   mkdir -p $OUTDIR

   foldseek easy-search $QUERYDIR $TARGETDIR $FSOUTFILE tmp --format-output "query,target,evalue,bits,lddt,prob,qcov,tcov,alntmscore" --exhaustive-search -e 10 --threads 1 -v 0

   OUTFILE="$OUTDIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10_mean.bitscore.tsv"
   python orthogroups/2_fully_merged_orthogroups/foldseek/parse_foldseek_query_target_pair.py $FSOUTFILE $OUTFILE $PDBDIR
   ```

   Prefilter to select candidates for the second merge.
   ```bash
   Rscript orthogroups/2_fully_merged_orthogroups/foldseek/get_hhsearch_candidates_for_foldseek_merge.R
   ```


   # Refined orthogroups from mapping singleton proteins, added proteins, and fusion proteins
   
   Map singleton proteins, added proteins, and fusion proteins to their fully merged orthogroups to generate refined orthogroups. This step takes compute-intensive hmmsearch results as inputs; precomputed inputs are available in ```data/orthogroups/```.
   ```bash
   Rscript orthogroups/3_refined_orthogroups/get_refined_orthogroups.R
   ```

   Resampled orthogroups by identifying subsets of protein members with conserved Pfam domains. This was used to resample two orthogroups in the study.
   ```bash
   Rscript orthogroups/3_refined_orthogroups/get_resampled_orthogroups.R
   ```

