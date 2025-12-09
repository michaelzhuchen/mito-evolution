#!/bin/bash

# Run DIAMOND commands from taskarray.

export IFS=$'\t'

for cmd in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
    outfile_prefix=`echo $cmd | cut -d" " -f 9` # extract the output filename from the -o argument of the command
    outfilename="${outfile_prefix}.gz" # add a .gz for compression
    
    # check to see if the file is already present. Only execute if it's absent.
    if [[ ! -s $outfilename ]]
    then
       eval "$cmd"; \
    fi
done;
