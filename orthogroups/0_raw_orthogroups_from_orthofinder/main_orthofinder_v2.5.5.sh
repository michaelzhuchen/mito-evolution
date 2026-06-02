#!/bin/bash

# Get raw orthogroups using precomputed Diamond results
conda activate OrthoFinder_v255
orthofinder -b $Orthofinder_Working_Directory_With_Diamond_Results -t 256 -a 32 -og