#!/bin/bash

EXEPATH=/usr/bin
EXELIST=( blastp blast_formatter $CAPS extractalign fastafetch fastaindex Gblocks mafft makeblastdb muscle parallel phyml prank seqkit squizz treebest )

check_bin(){
  for e in ${EXELIST[@]} ; do
    if [ -x "$EXEPATH/$e" ] ; then
      echo -e "[\e[92mOK\e[39m] $e"
    else
      echo -e "[\e[91mERROR\e[39m] $e not found!"
    fi
  done
}
