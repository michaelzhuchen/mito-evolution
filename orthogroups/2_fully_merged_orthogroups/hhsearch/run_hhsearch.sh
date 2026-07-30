#!/bin/sh

OGID=$1
QUERYFILE_SUFFIX=$2
QUERYDIR=$3
TARGETDB=$4
OUTDIR=$5

MSAFILE="${OGID}${QUERYFILE_SUFFIX}"

# Specify output file suffix
OUTFILENAME="${MSAFILE}_global.hhr"

OUTFILE="$OUTDIR/$OUTFILENAME"
if [[ ! -s "$OUTFILE" ]]
then
	# For global alignment, filter by probability>30 for hmm+structure merge
	hhsearch -i $QUERYDIR/$MSAFILE -d $TARGETDB -M 50 -maxres 65535 -p 30 -B 0 -b 0 -Z 500 -z 1 -o $OUTFILE -v 0 -glob
fi

# Process the HHR results to TSV
OUTFILE_TSV="${OUTDIR_TSV}/${OUTFILENAME}_prob30.tsv" # for hmm+structure merge

if [[ ! -s "$OUTFILE_TSV" ]]
then
	# Use python script to extract out the hits below the expect threshold
	python orthogroups/2_fully_merged_orthogroups/hhsearch/process_hhr_to_tsv_expect_threshold.py $OUTFILE $OUTFILE_TSV

	# Add the query OG as the first column
	sed -i "s/^/$OGID\t/" $OUTFILE_TSV
fi

