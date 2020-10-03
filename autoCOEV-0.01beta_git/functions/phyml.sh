#!/bin/bash

TREES="Trees"

# Construct the folder where trees will be calculated
MSATRETAG="${MSAMETHOD}_${PHYMLGBLOCKS}-${PHYMLGUIDE}-${TREESROOT}"
PHYMLCORES="--j ${THREADS}" # "12"

tree_names(){
  mkdir -p $TMP/$TREES
  cp $CWD/$EXTTREE $TMP/$TREES
  cat $CWD/$SPECIES | while read -r a b ; do
    sed -i "s/$b/$a/g" $TMP/$TREES/$EXTTREE
  done
  echo -e "Species names in the external tree changed to TAXID!"
}

# Trim the external tree following the species in each MSA. These can be
# used as guide trees for PhyML or as 'external' trees in the CAPS run.
ext_trees(){
  mkdir -p $TMP/$TREES/external
  cd $TMP/$MSA/$MSAMETHOD
  VANILLA=$( ls *.fa )
  for van in ${VANILLA[@]} ; do
    seqkit --threads ${THREADS} seq $van -n > ${van%.*}.species
  done
  SUBSPEC=$( ls *.species )
  for spcs in ${SUBSPEC[@]} ; do
    treebest subtree $TMP/$TREES/$EXTTREE $spcs > $TMP/$TREES/external/${spcs%.*}.tre
    echo -e "Trim external tree for $spcs"
  done
}

# Copy the necessary MSA files, converting the into PHYLIPI format.
phyml_prep() {
  mkdir -p $TMP/$TREES/phyml/$MSATRETAG/{phy,tre}
  cd $TMP/$MSA/$MSAMETHOD
    if [ "$GBLOCKS" = "vanilla" ]; then
      VANMSA=$(ls *.fa)
      for vmsa in ${VANMSA[@]} ; do
        squizz $vmsa -c PHYLIPI > $TMP/$TREES/phyml/$MSATRETAG/phy/${vmsa%.*}.phy
      done
    elif [ "$GBLOCKS" = "gblocks" ]; then
      GBLMSA=$(ls *.gbl)
      for gmsa in ${GBLMSA[@]} ; do
        squizz $gmsa -c PHYLIPI > $TMP/$TREES/phyml/$MSATRETAG/phy/${gmsa%.*.*}.phy
      done
    else
      echo -e "Check your MSAs!"
    fi
  cd $TMP/$TREES/phyml/$MSATRETAG/phy/
  ORTHPHY=$( ls *.phy )
  cd $TMP/$TREES/phyml/$MSATRETAG/tre
}

phymlfn() {
  local PTRE="${1}"
  
  if [ "$PHYMLGUIDE" = "exguide" ]; then
    local ETRE="--inputtree $TMP/$TREES/external/${PTRE%.*}.tre"
    echo -e "PhyML will use $TMP/$TREES/external/${PTRE%.*}.tre!"
  elif [ "$PHYMLGUIDE" = "noguide" ]; then
    ETRE=""
    echo -e "PhyML will not use a guide tree!"
  else
    echo -e "Check you input tree settings!"
  fi

  phyml -d aa \
    $PHYMLOPTIONS $ETRE \
    --no_memory_check \
    --leave_duplicates \
    -i $TMP/$TREES/phyml/$MSATRETAG/phy/$PTRE
}

phyml_results() {
  cd $TMP/$TREES/phyml/$MSATRETAG/phy/
  TREESOUT=$( ls *.phy_phyml_tree.txt )
  for t in ${TREESOUT[@]} ; do
    tname=$(basename $t .phy_phyml_tree.txt)
    
    if [ "$TREESROOT" = "rooted" ]; then
      treebest root $tname.phy_phyml_tree.txt > $TMP/$TREES/phyml/$MSATRETAG/tre/$tname.tre
    elif [ "$TREESROOT" = "noroot" ]; then
      cp $tname.phy_phyml_tree.txt $TMP/$TREES/phyml/$MSATRETAG/tre/$tname.tre
    else
      echo -e "\nCheck your settings!"
    fi
  
  done
}
