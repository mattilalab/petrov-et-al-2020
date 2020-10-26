#!/bin/bash

MSA="MSA" # MSA dir

CORESPRANK="--j ${THREADS}"
CORESMUSCLE="--j ${THREADS}"
THREADMAFFT="--thread ${THREADS}"

musclefn() {
  local SEQM="${1}"
  $MSAMETHOD $MUSCLEOPTIONS > $TMP/$MSA/$MSAMETHOD/${SEQM} < $SEQM
}

prankfn() {
  if [ "$PRANKGUIDE" = "exguide" ]; then
    mkdir -p $TMP/$MSA
    cp $CWD/$EXTTREE $TMP/$MSA
    cat $CWD/$SPECIES | while read -r a b ; do
    sed -i "s/$b/$a/g" $TMP/$MSA/$EXTTREE
    done
    echo -e "Species names in the external tree changed to TAXID!"
    PRGUID="-t=${TMP}/$MSA/${EXTTREE} -prunetree"
  else
    PRGUID=""
  fi
  
  local SEQP="${1}"
  $MSAMETHOD $PRGUID $PRANKOPTIONS \
    -o=$TMP/$MSA/$MSAMETHOD/$SEQP \
    -d=$SEQP
  mv $TMP/$MSA/$MSAMETHOD/$SEQP.best.fas $TMP/$MSA/$MSAMETHOD/$SEQP
}

# MAFFT has a nice multithreading option, which we take advantage of
mafftfn() {
  for seq in * ; do
    $MSAMETHOD $THREADMAFFT $MAFFTOPTIONS $seq > $TMP/$MSA/$MSAMETHOD/$seq
  done
}

# Run Gblocks and create a list of species of each MSA. Some help from:
# https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
# Also, make a list of species for each MSA
msa_process(){
  cd $TMP/$MSA/$MSAMETHOD
  PROCESSMSA=$(ls *.fa)
  for van in ${PROCESSMSA[@]} ; do
    Gblocks $van -t=p $GBLOCKSOPT -e=".gbl"
    sed -i "s: ::g" $van.gbl
    seqkit --threads ${THREADS} seq $van -n > ${van%.*}.species
  done
}
