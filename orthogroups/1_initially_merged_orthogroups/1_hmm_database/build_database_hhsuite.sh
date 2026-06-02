#!/bin/bash

BASEDIR="first_merge"
DBNAME="OG_raw_trimmed_hhsuite"

OUTFOLDER="${BASEDIR}/${DBNAME}_db"
N_CORES="8"

mkdir -p $OUTFOLDER

cd $OUTFOLDER

# Make index for MSAs
ffindex_build -s ${DBNAME}_msa.ff{data,index} $BASEDIR/$DBNAME

# Convert to a3m format
OMP_NUM_THREADS=${N_CORES} mpirun -np ${N_CORES} ffindex_apply_mpi ${DBNAME}_msa.ff{data,index} \
  -i ${DBNAME}_a3m_wo_ss.ffindex -d ${DBNAME}_a3m_wo_ss.ffdata \
    -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0
rm ${DBNAME}_msa.ff{data,index}

# Skip secondary structure prediction by renaming
mv ${DBNAME}_a3m_wo_ss.ffdata ${DBNAME}_a3m.ffdata
mv ${DBNAME}_a3m_wo_ss.ffindex ${DBNAME}_a3m.ffindex

# Build hmms
mpirun -np ${N_CORES} ffindex_apply_mpi ${DBNAME}_a3m.ff{data,index} -i ${DBNAME}_hhm.ffindex -d ${DBNAME}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0

# Computing context states for prefiltering
mpirun -np ${N_CORES} cstranslate_mpi -f -x 0.3 -c 4 -I a3m -i ${DBNAME}_a3m -o ${DBNAME}_cs219 

# Sort according to the number of columns in the MSAs
sort -k3 -n -r ${DBNAME}_cs219.ffindex | cut -f1 > sorting.dat
    
ffindex_order sorting.dat ${DBNAME}_hhm.ff{data,index} ${DBNAME}_hhm_ordered.ff{data,index}
mv ${DBNAME}_hhm_ordered.ffindex ${DBNAME}_hhm.ffindex
mv ${DBNAME}_hhm_ordered.ffdata ${DBNAME}_hhm.ffdata
    
ffindex_order sorting.dat ${DBNAME}_a3m.ff{data,index} ${DBNAME}_a3m_ordered.ff{data,index}
mv ${DBNAME}_a3m_ordered.ffindex ${DBNAME}_a3m.ffindex
mv ${DBNAME}_a3m_ordered.ffdata ${DBNAME}_a3m.ffdata
