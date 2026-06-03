#!/bin/sh

QUERY_FASTA="data/deeploc/swissprot_mitoepi_combined.fasta"
TARGET_FASTA="data/deeploc/swissprot_mitoepi_combined.fasta"
OUTPUT_FILE="data/deeploc/ggsearch/ggsearch_swissprot_mitoepi_combined_fasta_vs_self_expect1_combined.out"

fasta-36.3.8i/bin/ggsearch36 -E 1 -m 8 -T 2 $QUERY_FASTA $TARGET_FASTA > $OUTPUT_FILE