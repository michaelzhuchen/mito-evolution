#!/bin/bash

OUTPUT_DIRECTORY="output_directory"
TEMPORARY_DIRECTORY="tmp"

OGID="MOG0001047"

# Run MMseqs2
mmseqs easy-cluster ${OGID}.faa ${OUTPUT_DIRECTORY} ${TEMPORARY_DIRECTORY} --threads 1 --split-memory-limit 2G --min-seq-id 0.5 --cov-mode 2 -c 0.8 -v 1
