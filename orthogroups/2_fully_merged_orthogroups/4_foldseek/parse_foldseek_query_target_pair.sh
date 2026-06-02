#!/bin/sh

DATADIR="output_directory"
OUTDIR="output_directory"

QUERY_OG="MOG0030340"
TARGET_OG="OG0059097"


INPUTFILE="$DATADIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10.tsv"
OUTFILE="$OUTDIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10_mean.bitscore.tsv"

if [[ -s "$INPUTFILE" && ! -f "$OUTFILE" ]]
then
    python parse_foldseek_query_target_pair.py $INPUTFILE $OUTFILE $PDBDIR
fi