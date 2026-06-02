#!/bin/sh

OGID="MOG0030340"
OUTDIR_TSV="output_directory"
DB_NAME="second_merge/OG_merged_trimmed_hhsuite_db"
DBNAME="OG_merged_trimmed_hhsuite"


FAAFILE="$OGID"
MSAFILE="${FAAFILE}${QUERYFILE_SUFFIX}"

# Specify output file suffix
OUTFILENAME="${MSAFILE}_global.hhr"

OUTFILE="$OUTDIR/$OUTFILENAME"
if [[ ! -s "$OUTFILE" ]]
then
	# For global alignment, filter by probability>30 for hmm+structure merge
	hhsearch -i $FAADIR/$MSAFILE -d $DB_NAME -M 50 -maxres 65535 -p 30 -B 0 -b 0 -Z 500 -z 1 -o $OUTFILE -v 0 -glob
fi


# Process the HHR results to TSV
OUTFILE_TSV="${OUTDIR_TSV}/${OUTFILENAME}_prob30.tsv" # for hmm+structure merge

if [[ ! -s "$OUTFILE_TSV" ]]
then
	# Use python script to extract out the hits below the expect threshold
	python process_hhr_to_tsv_expect_threshold.py $OUTFILE $OUTFILE_TSV

	# Add the query OG as the first column
	sed -i "s/^/$OGID\t/" $OUTFILE_TSV
fi

