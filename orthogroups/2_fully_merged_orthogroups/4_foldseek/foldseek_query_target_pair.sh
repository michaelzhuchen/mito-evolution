#!/bin/sh

DATADIR="output_directory"
OUTDIR="output_directory"
TEMPDIR="tmp"

QUERY_OG="MOG0030340"
TARGET_OG="OG0059097"


OUTFILE="$OUTDIR/${QUERY_OG}_vs_${TARGET_OG}_exhaustive_expect10.tsv"

if [[ ! -f $OUTFILE ]]
then
    QUERYFILE="$DATADIR/${QUERY_OG}/structs_raw"
    TARGETFILE="$DATADIR/${TARGET_OG}/structs_raw"
    TEMPFOLDER="$TEMPDIR/${QUERY_OG}_${TARGET_OG}"
    mkdir -p $TEMPFOLDER

    foldseek easy-search $QUERYFILE $TARGETFILE $OUTFILE $TEMPFOLDER --format-output "query,target,evalue,bits,lddt,prob,qcov,tcov,alntmscore" --exhaustive-search -e 10 --threads 1 -v 0
fi
