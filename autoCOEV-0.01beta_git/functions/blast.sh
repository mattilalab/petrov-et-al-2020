#!/bin/bash

GETFA="BestBLASTfasta" # Best pBLAST hits

# Retrieve the mouse protein sequences for each gene from UniProt and store them separately.
# Retry up to 50 times (or set "-t inf" for infinite retries). Check for empty files if
# download was unsuccessful. I can parallelize this, but the server may block the download.
uniprot_download() {
  cd $TMP/$ORTHO/
  for UniProt in * ; do
  mkdir -p $UniProt/$ORGANISM
  wget \
    -t 50 \
    -c "https://www.uniprot.org/uniprot/${UniProt}.fasta" \
    -O $UniProt/$ORGANISM/$UniProt.fa
  done
  find $TMP/$ORTHO/*/$ORGANISM/ -size 0 -print >> $TMP/tsv/UniProt.failed
}

# Create BLAST database for each mouse sequence, against which we will run the sequences of
# collected orthologues in order to determine best matching isoform. Use only the found proteins.
blast_db_prep() {
  local UniProt=${1}
  makeblastdb \
    -in $UniProt/$ORGANISM/$UniProt.fa \
    -parse_seqids \
    -blastdb_version 5 \
    -title "$UniProt" \
    -dbtype prot
}

# Run reciprocal BLASTP against mouse protein and sort by bit score. Then sort by "bit score"
# (12th column, --key 12). Output format info: http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
reciprocal_blast() {
  local UniProt=${1}
  cd $UniProt/FASTA/
  for blfas in * ; do
    mkdir -p ../BLAST/${blfas%.*}
    blastp \
      -query $blfas \
      -db ../$ORGANISM/$UniProt.fa \
      -out ../BLAST/${blfas%.*}/${blfas%.*}.out \
      -outfmt="6 qseqid sacc pident ppos length mismatch gapopen gaps qstart qend sstart send evalue bitscore"
    sort --key 14 --numeric-sort --reverse ../BLAST/${blfas%.*}/${blfas%.*}.out | head -n 1 >> ../blastBest.tsv
    echo -e "BLAST ${blfas%.*} against $UniProt"
    
    # Shall we run BLAST again to generate detailed results?
    if [ "$DETBLAST" = "yes" ]; then
      blastp \
        -query $blfas \
        -db ../$ORGANISM/$UniProt.fa \
        -out ../BLAST/${blfas%.*}/${blfas%.*}.default \
        -outfmt="0"
     echo -e "detailed..."
     elif [ "$DETBLAST" = "no" ]; then
       return 1
    else
      "Check your BLAST settings!"
    fi
  done
  cd ../../
  echo ""
}

# Get the fasta sequences of the best hits with certain identity, defined by user. Help from:
# https://stackoverflow.com/questions/8654051/how-to-compare-two-floating-point-numbers-in-bash
best_hits() {
  local UniProt=${1}
  mkdir -p $TMP/$GETFA
  cat $UniProt/blastBest.tsv | \
  while read -r qseqid sacc pident ppos length mismatch gapopen gaps qstart qend sstart send evalue bitscore ; do
    PERCENTGAPS=$(echo "${gaps}/${length}*100" |bc -l)
    if (( $(echo "$pident > $PIDENT" |bc -l) && $(echo "$PERCENTGAPS < $PGAPS" |bc -l) )); then
      echo -e "[\e[92mGOOD\e[39m] Identity and gaps: $sacc <- $qseqid: $pident and $PERCENTGAPS"
      fastafetch \
        -f $DTB/$ALLFASTA \
        -i $DTB/$ALLFASTA.index \
        -q $qseqid >> $TMP/$GETFA/$UniProt.fa
      echo -e "$qseqid\t$sacc\t$pident\t$ppos\t$length\t$mismatch\t$gapopen\t$gaps\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore\t$PERCENTGAPS" >> $TMP/tsv/blastBestFasta.tsv
    else
      echo -e "[\e[91mPOOR\e[39m] Identity and gaps: $sacc <- $qseqid: $pident and $PERCENTGAPS"
      echo -e "$qseqid\t$sacc\t$pident\t$ppos\t$length\t$mismatch\t$gapopen\t$gaps\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore\t$PERCENTGAPS" >> $TMP/tsv/blastBestExclude.tsv
    fi
  done
}

# Replace fasta titles with species latin names. In the end, add the
# sequence from UniProt
species_names() {
  local UniProt=${1}
  rootUP=$( basename $UniProt .fa )
  sed -i "s/_.*//" $UniProt
  echo ">${ORGANISM}" >> $UniProt
  sed 1d $TMP/$ORTHO/$rootUP/$ORGANISM/$UniProt >> $UniProt
}

headers_blast(){
  while read -r UniProtID geneName assign ; do
    sed -i "s/$UniProtID/$UniProtID\t$geneName/g" $TMP/tsv/blastBestFasta.tsv
    sed -i "s/$UniProtID/$UniProtID\t$geneName/g" $TMP/tsv/blastBestExclude.tsv
  done < $CWD/$PROTEIN
  
  while read -r taxidID speciesName ; do
    sed -i "s/${taxidID}_/$speciesName\t${taxidID}\t${taxidID}_/g" $TMP/tsv/blastBestFasta.tsv
    sed -i "s/${taxidID}_/$speciesName\t${taxidID}\t${taxidID}_/g" $TMP/tsv/blastBestExclude.tsv
  done < $CWD/$SPECIES
  
  sed -i "1i Species\ttaxid\tqseqid\tsacc\tname\tpident\tppos\tlength\tmismatch\tgapopen\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tgaps" $TMP/tsv/blastBestFasta.tsv
  sed -i "1i Species\ttaxid\tqseqid\tsacc\tname\tpident\tppos\tlength\tmismatch\tgapopen\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tgaps" $TMP/tsv/blastBestExclude.tsv
}
