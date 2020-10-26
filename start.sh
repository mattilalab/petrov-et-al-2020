#!/bin/bash

# Export all
set -a

# What date and time is it?
DATESTAMP=$(date)

# Current work dir
CWD=$(pwd)

NAME="AutoCoEv" # script name
VER="0.021beta" # version

# Load settings
. $CWD/settings.conf

# Load functions
. $CWD/functions/databases.sh
. $CWD/functions/retrieval.sh
. $CWD/functions/blast.sh
. $CWD/functions/msa.sh
. $CWD/functions/trees.sh
. $CWD/functions/pairing.sh
. $CWD/functions/caps.sh
. $CWD/functions/results.sh
. $CWD/functions/xml.sh
. $CWD/functions/check.sh

# User provided input
echo -e "\n\e[46mWelcome to $NAME version $VER!\e[49m\n"

# Check for binaries
echo -e "Checking for executables in: \e[96m$EXEPATH\e[39m"
check_bin
echo ""

read -p "Change working dir or press ENTER to accept: " -e -i $TMP TMP
echo -e "\n\e[96m$TMP\e[39m\n"

echo -e "Do initial preparations:"
# Help from:
# https://unix.stackexchange.com/questions/293340/bash-how-can-i-re-display-selection-menu-after-a-selection-is-chosen-and-perfo
PS3="Your choice: "
PREPOPT=( "Download databases"
         "Extract databases"
         "Index FASTA database"
         "Trim gene_xrefs database ($ORGANISM)"
	 "Trim OG2genes database ($LEVEL)"
         "[DONE AND CONTINUE]" )

select prep in "${PREPOPT[@]}" ; do
case $prep in

"Download databases")
  download_db
echo -e "\nDone with 1)"
;;

"Extract databases")
  extract_db
echo -e "\nDone with 2)"
;;

"Index FASTA database")
  index_fa
echo -e "\nDone with 3)"
;;

"Trim gene_xrefs database ($ORGANISM)")
  trim_gene_xrefs
echo -e "\nDone with 4)"
;;

"Trim OG2genes database ($LEVEL)")
  trim_gene_xrefs
echo -e "\nDone with 5)"
;;

"[DONE AND CONTINUE]")
echo -e "Continue...\n"
break 2
;;

*)
echo "Invalid option"
;;

esac
REPLY=
done

cd $CWD

read -p "Gene XRefs database: " -e -i $GENEXREF GENEXREF
echo -e "\n\e[96m$GENEXREF\e[39m"
head $DTB/$GENEXREF
echo ""

read -p "OrthoGroup 2 genes database: " -e -i $OG2GENES G2GENES
echo -e "\n\e[96m$OG2GENES\e[39m"
head $DTB/$OG2GENES
echo ""

read -p "All FASTA sequences database: " -e -i $ALLFASTA ALLFASTA
echo -e "\n\e[96m$ALLFASTA\e[39m"
seqkit head $DTB/$ALLFASTA
echo -e "\n\e[96m$ALLFASTA.index\e[39m"
head $DTB/$ALLFASTA.index
echo ""

read -p "Change proteins list or press ENTER to accept: " -e -i $PROTEIN PROTEIN
echo -e "\n\e[96m$PROTEIN\e[39m"
head $PROTEIN
echo ""

read -p "Change species list or press ENTER to continue: " -e -i $SPECIES SPECIES
echo -e "\n\e[96m$SPECIES\e[39m"
head $SPECIES
echo ""

read -p "Use maximum available CPU threads?: " -e -i $THREADS THREADS
echo -e "CPU threads to be used are: \e[92m$THREADS\e[39m\n"
echo -e "We achieve parallelization via \e[92mGNU/Parallel\e[39m"
echo -e "Run 'parallel --citation' once to silence its citation notice.\n\n"

# Output a summary of settings
echo -e "Work directory: \e[92m${TMP}\e[39m\n"
echo -e "UniProt reference sequences are from:.........\e[92m$ORGANISM\e[39m"
echo -e "Reverse BLAST identity (%) cutoff:............\e[92m$PIDENT\e[39m"
echo -e "Reverse BLAST gaps (%) cutoff:................\e[92m$PGAPS\e[39m"
echo -e "MSAs for PhyML and/or CAPS created by:........\e[92m$MSAMETHOD\e[39m"
echo -e "Process MSAs before CAPS run (Gblocks)?.......\e[92m$GBLOCKS\e[39m"
echo -e "Run PhyML on processed MSAs (Gblocks)?........\e[92m$PHYMLGBLOCKS\e[39m"
echo -e "PhyML user-specified options:.................\e[92m$PHYMLOPTIONS\e[39m"
echo -e "Root the PhyML produced trees?................\e[92m$TREESROOT\e[39m"
echo -e "Do we use external guide tree for PhyML?......\e[92m$PHYMLGUIDE\e[39m"
echo -e "Trees to be used with CAPS:...................\e[92m$TREESCAPS\e[39m"
echo -e "Minimum number of common species in a pair:...\e[92m$MINCOMMONSPCS\e[39m"
echo -e "CAPS alpha-value cutoff at runtime:...........\e[92m$ALPHA\e[39m"
echo -e "CAPS bootstrap value at runtime:..............\e[92m$BOOT\e[39m"
echo -e "Postrun bootstrap threshold cutoff:...........\e[92m$RESBOOT\e[39m"
echo -e "Alignment gaps postrun cutoff:................\e[92m$RESGAPS\e[39m"
echo -e "Column identity postrun cutoff:...............\e[92m$RESIDEN\e[39m"
echo -e "Postrun P-value correlation cutoff:...........\e[92m$PVALUE\e[39m"
echo -e "\n"

echo -e "Select a step:"

mkdir -p $TMP/tsv
cd $TMP

# Menu. Great help from:
# https://askubuntu.com/questions/1705/how-can-i-create-a-select-menu-in-a-shell-script

PS3="Your choice: "
SEQOPT=( "Pair UniProt <-> OrthoDB <-> OGuniqueID"
         "Prepare orthologues list (level: $LEVEL)"
         "Get FASTA sequences of all orthologues"
         "Download sequences from UniProt (organism: $ORGANISM)"
         "BLAST orthologues against UniProt sequence ($ORGANISM, detailed: $DETBLAST)"
         "Get FASTA sequences of the best hits (identity: $PIDENT; gaps: $PGAPS)"
         "[MSA] Create MSA with selected method ($MSAMETHOD)"
         "[TRE] Prepare trees ($TREESCAPS, $MSAMETHOD, $PHYMLGBLOCKS, $PHYMLGUIDE, $TREESROOT)"
         "[RUN] Create pairs ($PAIRINGMANNER)"
         "[RUN] CAPS run (alpha: $ALPHA, $MSAMETHOD, $GBLOCKS, $TREESCAPS)"
         "[RES] Inspect CAPS results"
         "[RES] Generate columns stats"
         "[XML] Process CAPS results"
         "[Exit script]" )

select opt in "${SEQOPT[@]}" ; do
case $opt in

"Pair UniProt <-> OrthoDB <-> OGuniqueID")
  create_tsv_files
  pair_uniprot_vs_orthodb

  cd $TMP/$ORTHO/
  FOUND=$( ls ./ )
  parallel $CORESCAPS extract_oguniqueid ::: $FOUND
  parallel $CORESCAPS report_gen ::: $FOUND
  cd ..
  
  summary_clarify
  report_duplicates
  headers_insert
echo -e "\nDone with 1)"
;;

"Prepare orthologues list (level: $LEVEL)")
  sed "/$ORGANISM/d" $CWD/$SPECIES | \
  while read -r taxidNCBI latinName; do
    cd $TMP/$ORTHO/
    FOUND=$( ls ./ )
    parallel $CORESCAPS prepare_orthologues_list ::: $FOUND
    cd ..
  done
  
  cd $TMP/$ORTHO/
  FOUND=$( ls ./ )
  parallel $CORESCAPS no_homologues_check ::: $FOUND
  cd ..
echo -e "\nDone with 2)"
;;

"Get FASTA sequences of all orthologues")
  cd $TMP/$ORTHO/
  FOUND=$( ls ./ )
  parallel $CORESCAPS get_ortho_fasta ::: $FOUND
  cd ..
echo -e "\nDone with 3)"
;;

"Download sequences from UniProt (organism: $ORGANISM)")
  uniprot_download
echo -e "\nDone with 4)"
;;

"BLAST orthologues against UniProt sequence ($ORGANISM, detailed: $DETBLAST)")
  cd $TMP/$ORTHO/
    FOUND=$( ls ./ )
    parallel $CORESCAPS blast_db_prep ::: $FOUND
    parallel $CORESCAPS reciprocal_blast ::: $FOUND
  cd ..
echo -e "\nDone with 5)"
;;

"Get FASTA sequences of the best hits (identity: $PIDENT; gaps: $PGAPS)")
  cd $TMP/$ORTHO/
  FOUND=$( ls ./ )
  parallel $CORESCAPS best_hits ::: $FOUND
  cd ..

  cd $TMP/$GETFA/
  FOUND=$( ls ./ )
  parallel $CORESCAPS species_names ::: $FOUND
  cd ..
  
  headers_blast
  echo -e "\nDone with 6)"
;;

"[MSA] Create MSA with selected method ($MSAMETHOD)")
  mkdir -p $TMP/$MSA/$MSAMETHOD/
  cd $TMP/$GETFA
  ORTHMSA=$( ls ./ )
  
  if [ "$MSAMETHOD" = "muscle" ]; then
    parallel $CORESMUSCLE musclefn ::: "$ORTHMSA"
  elif [ "$MSAMETHOD" = "prank" ]; then
    parallel $CORESPRANK prankfn ::: "$ORTHMSA"
  elif [ "$MSAMETHOD" = "mafft" ] || [ "$MSAMETHOD" = "mafft-linsi" ] || [ "$MSAMETHOD" = "mafft-ginsi" ] || [ "$MSAMETHOD" = "mafft-einsi" ] || [ "$MSAMETHOD" = "mafft-fftns" ] || [ "$MSAMETHOD" = "mafft-fftnsi" ]; then
    mafftfn
  else
    echo "[ERROR] Specify MSA method properly!"
  fi
  cd ..
  msa_process
echo -e "\nDone with 7)"
;;

"[TRE] Prepare trees ($TREESCAPS, $MSAMETHOD, $PHYMLGBLOCKS, $PHYMLGUIDE, $TREESROOT)")
if [ "$TREESCAPS" = "auto" ]; then
  echo -e "CAPS will generate its own trees. Skipping..."
elif [ "$TREESCAPS" = "external" ]; then
  ext_trees
elif [ "$TREESCAPS" = "phyml" ] && [ "$PHYMLGUIDE" = "exguide" ]; then
  phyml_guide
  phyml_prep
  parallel $CORESPHYML phymlfn ::: "$ORTHPHY"
  phyml_process
elif [ "$TREESCAPS" = "phyml" ] && [ "$PHYMLGUIDE" = "noguide" ]; then
  phyml_prep
  parallel $CORESPHYML phymlfn ::: "$ORTHPHY"
  phyml_process
else
  echo -e "Check your tree settings!"
fi
echo -e "\nDone with 8)"
;;

"[RUN] Create pairs ($PAIRINGMANNER)")
if [ "$PAIRINGMANNER" = "all" ]; then
  pair_msa
  pair_tree
elif [ "$PAIRINGMANNER" = "defined" ]; then
  pair_defined_msa
  pair_tree
else
  echo -e "Check your pairing settings!"
fi
echo -e "\nDone with 9"
;;

"[RUN] CAPS run (alpha: $ALPHA, $MSAMETHOD, $GBLOCKS, $TREESCAPS)")
split_dirs
caps_set
for folder in $TMP/$CAPSM/* ; do
  cd $folder
  echo -e "Processing \e[92m${folder}\e[39m"
  echo "$DATESTAMP ${folder}" >> $TMP/progress-$ALPHA-$MSAMETHOD-$GBLOCKS-$TREESCAPS.txt
  PAIRLIST=$(ls)
  parallel $CORESCAPS --progress capsfn ::: "$PAIRLIST"
  echo -e "Done in \e[92m${folder}\e[39m"
  echo "" >> $TMP/progress-$ALPHA-$MSAMETHOD-$GBLOCKS-$TREESCAPS.txt
  cd ..
done
echo -e "\nDone with 10"
;;

"[RES] Inspect CAPS results")
cd $TMP/$CAPSM/
for folder in * ; do
  cd $folder
  SUBPART=$( ls ./ )
  parallel $CORESCAPS caps_inspect ::: "$SUBPART"
  cd ..
done
echo -e "\n\e[92mResults inspections done!\e[39m\n"

cd $TMP/$RESULTS/coev/
for folder in * ; do
  cd $folder
  SUBPART=$( ls ./ )
  parallel $CORESCAPS results_cleanup ::: "$SUBPART"
  cd ..
done

echo -e "\nDone with 11"
;;

"[RES] Generate columns stats")

cd $TMP/$RESULTS/coev/
for folder in * ; do
  cd $folder
  SUBPART=$( ls ./ )
  parallel $CORESCAPS columns_stats ::: "$SUBPART"
  cd ..
done

cd $TMP/$RESULTS/coev/
for folder in * ; do
  cd $folder
  SUBPART=$( ls ./ )
  parallel $CORESCAPS mean_pval ::: "$SUBPART"
  cd ..
done

summary_cleanup

echo -e "\nDone with 12"
;;

"[XML] Process CAPS results")
node_gen_list
make_coev_inter_xml
make_pos_corr_xml
echo -e "\nDone with 13"
;;

"[Exit script]")
echo "Quitting..."
exit
;;

*)
echo "Invalid option"
;;

esac
REPLY=
done
