#!/bin/bash

ALLFASTA="odb10v1_all_fasta.tab" # FASTA db
GENEXREF="odb10v1_gene_xrefs.${ORGANISM}.tab" # Genes from a certain organism (e.g. ORGANISM="10090" for mouse)
OG2GENES="odb10v1_OG2genes.${LEVEL}.tab" # orthologues groups db at specified level (e.g. LEVEL="32523")

download_db(){
  mkdir -p $DTB
  cd $DTB
  wget -c "https://v101.orthodb.org/download/odb10v1_all_fasta.tab.gz"
  wget -c "https://v101.orthodb.org/download/odb10v1_gene_xrefs.tab.gz"
  wget -c "https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz"
  echo -e "Downloaded the three essential databases in \e[96m$DTB\e[39m"
}

extract_db(){
  cd $DTB
  gunzip -v -c odb10v1_all_fasta.tab.gz > odb10v1_all_fasta.tab
  gunzip -v -c odb10v1_gene_xrefs.tab.gz > odb10v1_gene_xrefs.tab
  gunzip -v -c odb10v1_OG2genes.tab.gz > odb10v1_OG2genes.tab
  echo -e "Extracted the databases in \e[96m$DTB\e[39m. You may delete the archives."
}

index_fa(){
  cd $DTB
  fastaindex \
    odb10v1_all_fasta.tab \
    odb10v1_all_fasta.tab.index
  echo -e "Indexing of odb10v1_all_fasta.tab in \e[96m$DTB\e[39m complete!"
}

trimmed_db(){
  cd $DTB
  grep "\b${ORGANISM}_" odb10v1_gene_xrefs.tab > odb10v1_gene_xrefs.${ORGANISM}.tab
  grep "at${LEVEL}\b" odb10v1_OG2genes.tab > odb10v1_OG2genes.${LEVEL}.tab
  echo -e "Sub-databases extracted in \e[96m$DTB\e[39m . You may delete the original."
}
