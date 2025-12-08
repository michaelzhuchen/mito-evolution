#!/bin/sh

# run just the top part of the for loop to do only unidirectional search
for OGID in `grep -P -w "$SGE_TASK_ID$" $SH_TASKARRAY_FILE | cut -f 1`; do \
  mkdir -p $TEMPOUTDIR/$OGID
  
  original_msa_file="$BASEDIR/$OGID.faa_clipkit.gappy.msa"

  # Convert lta2019 ids to ltaref
  msa_file="$TEMPOUTDIR/$OGID/$OGID.faa_clipkit.gappy.msa"
  Rscript convert_msa_protein_ids.R $original_msa_file $msa_file

  # Use protein ids from PhROGs
  Rscript get_PhROG_protein_ids.R $PHROGS_DIR $TEMPOUTDIR/$OGID $OGID
  
  for protein_id_list in $TEMPOUTDIR/$OGID/*_LOO_protein_id_list.txt
  do
    PHROGID=${protein_id_list%_LOO_protein_id_list.txt}
    PHROGID=${PHROGID##*/}

    phrog_msa_file="$TEMPOUTDIR/$OGID/${PHROGID}.faa_clipkit.gappy.msa"

    # Slice MSA to get just the PhROG proteins
    awk -v protein_id_list="$protein_id_list" '
    BEGIN {
        while ((getline id < protein_id_list) > 0) {
            gsub(/^[ \t]+|[ \t]+$/, "", id)  # Trim spaces
            ids[id] = 1
        }
        close(protein_id_list)
    }

    /^>/ {
        split(substr($0, 2), parts, " ")
        seq_id = parts[1]
        keep = (seq_id in ids)
    }

    {
        if (keep) print
    }
    ' $msa_file > $phrog_msa_file

    hmmbuild -o /dev/null -O "$TEMPOUTDIR/$OGID/$PHROGID.hmm.msa" "$TEMPOUTDIR/$OGID/$PHROGID.hmm" "$phrog_msa_file"
    
    while IFS= read -r protein_id
    do
      OUTPUTPREFIX=$(echo "$TEMPOUTDIR/$OGID/${PHROGID}_LOO_${protein_id}")

      # Build hmm
      PROFILEHMM="$OUTPUTPREFIX.hmm"
      if [ ! -s "$PROFILEHMM" ]
      then
        awk -v seq="$protein_id" '
        BEGIN {keep=0}
        /^>/ {exclude=($0 ~ ">"seq) ? 1 : 0}
        !exclude
        ' "$phrog_msa_file" | \
        awk -v seq="$protein_id" '
        BEGIN {exclude=0}
        /^>/ {exclude=($0 ~ ">"seq) ? 1 : 0}
        !exclude
        ' "$phrog_msa_file" | \
        hmmbuild --informat afa -n "${PHROGID}_LOO_${protein_id}" --cpu 1 --seed 42 -o /dev/null -O "$PROFILEHMM.msa" "$PROFILEHMM" -
      fi

      # Hmmsearch
      ORIGINAL_FAAFILE="$BASEDIR/$OGID.faa"
      FAAFILE="$TEMPOUTDIR/$OGID/$OGID.faa"
      Rscript convert_fasta_protein_ids.R $ORIGINAL_FAAFILE $FAAFILE
      OUTFILE="${OUTPUTPREFIX}_hmmsearch.out"
      if [ ! -s "$OUTFILE" ]
      then
        awk -v seq="$protein_id" '
        BEGIN {include=0}
        /^>/ {include=($0 ~ ">"seq) ? 1 : 0}
        include
        ' "$FAAFILE" | \
        hmmsearch -o /dev/null --noali -Z $DBSIZE -E $EXPECT_THRESHOLD --tformat fasta --cpu 1 --domtblout "$OUTFILE" "$PROFILEHMM" -

        # Parse output to extract just the rows containing LOO protein id (remove hmmer comment lines)
        awk -v p="$protein_id" '$0 ~ ("^" p)' $OUTFILE > $OUTFILE.tmp && mv $OUTFILE.tmp $OUTFILE
      fi

      # Remove temporary files
      rm $PROFILEHMM
      rm $PROFILEHMM.msa

    done < "$protein_id_list"

    # Combine output files per PhROG
    cat $TEMPOUTDIR/$OGID/$PHROGID*_hmmsearch.out > $OUTDIR/results/hmmsearch_out/${PHROGID}_LOO_hmmsearch.out

  done
  
  # Clean up temporary files
  rm -r $TEMPOUTDIR/$OGID
done;

