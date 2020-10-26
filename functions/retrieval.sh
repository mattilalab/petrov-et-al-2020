#!/bin/bash

ORTHO="Orthologues"

# Create blank files, so the script does not complain if some was empty
# (e.g. if there were no duplicates)
create_tsv_files(){
  touch $TMP/tsv/OrthoDB_Missing.tsv
  touch $TMP/tsv/proteinsFound.tsv
  touch $TMP/tsv/Summary.tsv
  touch $TMP/tsv/duplicates_OrthoDB.tsv
  touch $TMP/tsv/duplicates_UniProt.tsv
  touch $TMP/tsv/duplicates_OrthoGroup.tsv
}

# Pair UniProt ID to the corresponding OrthoDB ID. Help with defining word borders with grep:
# https://www.linuxquestions.org/questions/linux-newbie-8/how-to-grep-for-an-exact-word-4175439257/
pair_uniprot_vs_orthodb(){
  while read -r UniProtID geneName assign ; do
    mkdir -p $TMP/$ORTHO/$UniProtID
    if LANG=C grep -F -w "${UniProtID}" $DTB/$GENEXREF | grep "\b${ORGANISM}_" >> $TMP/$ORTHO/$UniProtID/UniProt_OrthoDB.tsv ; then
      echo -e "[\e[92mMATCHED\e[39m] $UniProtID\t$geneName\t$assign"
    else
      echo -e "[\e[91mMISSING\e[39m] $UniProtID\t$geneName\t$assign"
      echo -e "$UniProtID\t$geneName\t$assign" >> $TMP/tsv/OrthoDB_Missing.tsv
      rm -rf $TMP/$ORTHO/$UniProtID
    fi
  done < $CWD/$PROTEIN
  echo ""
}

# Extract OGuniqueID. Use grep with LANG=C and -F option to speed up the process (no strings accepted).
extract_oguniqueid() {
  local UniProt="${1}"
  while read -r OrthoDBgeneID UniProtID externlDBname ; do
    if LANG=C grep -F -w "$OrthoDBgeneID" $DTB/$OG2GENES | grep "at${LEVEL}\b" >> $UniProt/OrthoGroup.tsv ; then
      echo -e "[\e[92mEXTRACT\e[39m] $UniProtID\t$OrthoDBgeneID"
    else
      echo -e "[\e[91mMISSING\e[39m] $UniProtID\t$OrthoDBgeneID"
    fi
  done < $UniProt/UniProt_OrthoDB.tsv
}

# Generate report. Some help with awk from:
# https://unix.stackexchange.com/questions/222121/how-to-remove-a-column-or-multiple-columns-from-file-using-shell-command
report_gen(){
  local UniProt="${1}"
  cd $UniProt
  paste OrthoGroup.tsv UniProt_OrthoDB.tsv | awk '!($3=$5="")'>> $TMP/tsv/Summary.tsv
  cd ..
}

summary_clarify(){
while read -r UniProtID geneName assign ; do
  sed -i "s:$UniProtID:$UniProtID $geneName $assign:g" $TMP/tsv/Summary.tsv
  sed -i "s:  : :g" $TMP/tsv/Summary.tsv
  sed -i "s: :	:g" $TMP/tsv/Summary.tsv
done < $CWD/$PROTEIN
}

# Deal with duplicate entries
# https://askubuntu.com/questions/434545/identify-duplicate-lines-in-a-file-without-deleting-them
# https://unix.stackexchange.com/questions/224433/grep-first-column-uniq-values/224434
# https://superuser.com/questions/1092282/bash-sort-by-not-first-character
report_duplicates(){

  # Output duplicated Ortho Groups entries
  for dup in $(awk '{ print $1 }' $TMP/tsv/Summary.tsv | sort | uniq -d) ; do
    grep "$dup" $TMP/tsv/Summary.tsv >> $TMP/tsv/duplicates_OrthoGroup.tsv
  done
  
  # Output duplicated OrthoDB entries
  for dup in $(awk '{ print $2 }' $TMP/tsv/Summary.tsv | sort | uniq -d) ; do
    grep "$dup" $TMP/tsv/Summary.tsv >> $TMP/tsv/duplicates_OrthoDB.tsv
  done

  # Output duplicated UniProt entries
  for dup in $(awk '{ print $3 }' $TMP/tsv/Summary.tsv | sort |  uniq -d) ; do
    grep "$dup" $TMP/tsv/Summary.tsv >> $TMP/tsv/duplicates_UniProt.tsv
  done
  
  # Output found proteins
  sort --key 3 -u $TMP/tsv/Summary.tsv | sort --key 4 | awk '{ print $3, $4, $5 }' >> $TMP/tsv/proteinsFound.tsv
  
  echo -e "Reporting duplicates..."
}

# Output dir and column headers. Do this in the end, so headers do not mess with the while loops
headers_insert(){
  sed "1i UniProt\tName\tAssign" -i $TMP/tsv/OrthoDB_Missing.tsv
  sed "1i UniProt\tName\tAssign" -i $TMP/tsv/proteinsFound.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/Summary.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/duplicates_OrthoDB.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/duplicates_UniProt.tsv
  sed "1i OGuniqueID\tOrthoDBunique\tUniProt\tName\tAssign" -i $TMP/tsv/duplicates_OrthoGroup.tsv
}

# Prepare the list of orthologues, collect all identifiers for all isoforms for each protein.
# Skip reference organism (e.g. mouse), as we will get the sequences from UniProt anyway.
prepare_orthologues_list() {
  local UniProt="${1}" 
  while read -r OGuniqueID OrthoDBgeneID ; do
    if LANG=C grep -F -w "$OGuniqueID" $DTB/$OG2GENES | grep "\<${taxidNCBI}_\B" >> $UniProt/orthologues.tsv ; then
      echo -e "[\e[92mMATCHED\e[39m] $OGuniqueID\t$UniProt\t$latinName"
      echo -e "$taxidNCBI\t$latinName" >> $UniProt/${OrthoDBgeneID}.speciesFound.tsv
    else
      echo -e "[\e[91mMISSING\e[39m] $OGuniqueID\t$UniProt\t$latinName"
      echo -e "$taxidNCBI\t$latinName" >> $UniProt/${OrthoDBgeneID}.speciesMissing.tsv
    fi
  done < $UniProt/OrthoGroup.tsv
}

# Check if no homologues were found. These will be removed
no_homologues_check(){
  local UniProt="${1}"
  find $TMP/$ORTHO/$UniProt/orthologues.tsv -size 0 -print >> $TMP/tsv/homologues_None.tsv \
  -exec rm -rf $TMP/$ORTHO/$UniProt {} \;
}

# Get the fasta sequences of orthologues all isoforms. Remove any duplicates.
# https://stackoverflow.com/questions/2664740/extract-file-basename-without-path-and-extension-in-bash
get_ortho_fasta() {
  local UniProt="${1}"
  cat $UniProt/orthologues.tsv | sort | uniq | \
  while read -r CLIDatCLADE OrthoDBgeneID ; do
    mkdir -p $UniProt/FASTA
    TAXID=$(basename "$OrthoDBgeneID" | cut -d_ -f1)
    fastafetch -s \
      -f $DTB/$ALLFASTA \
      -i $DTB/$ALLFASTA.index \
      -q $OrthoDBgeneID >> $UniProt/FASTA/$TAXID.fa
    echo -e "$UniProt\t$OrthoDBgeneID\t$TAXID"
  done
}
