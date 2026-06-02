#!/bin/bash

OGID="MOG0030340"

FAAFILE="${OGID}.faa"
MSAFILE="${FAAFILE}.msa"

# Align
conda activate mafft
mafft --auto --quiet --thread 1 --anysymbol $FAAFILE > $MSAFILE

# Trim
conda activate clipkit_2.3.0
clipkit -q -m gappy "${MSAFILE}" -o "${FAAFILE}_clipkit.gappy.msa"
python utils/remove_gap_only_seqs.py "${FAAFILE}_clipkit.gappy.msa"
