#!/bin/bash

# Inspect CAPS results. Separate pair folders where CAPS run failed,
# skip pair folders where no coevolution was detected and copy those
# where co-evolving pairs were found. To do this, first check if row
# 2 in coev_inter.csv is empty (failed pairs), then check if its third
# column value (coevolving sites) is equal to zero (no coevolution) or
# larger than zero (coevolution detected).
caps_inspect() {
local pair="${1}"
  SUMMARY=$(sed -n '2p' $pair/coev_inter.csv)
  while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
    if [ -z "${SUMMARY}" ]; then
      echo -e "[\e[91mFAILED\e[39m] Copying pair where CAPS failed: $pair"
      mkdir -p $TMP/$RESULTS/fail
      cp -a $pair $TMP/$RESULTS/fail
    elif [ "$numPairs" -eq 0 ]; then
      echo -e "[\e[34mNOCOEV\e[39m] No co-evolving pairs found for: $pair"
      mkdir -p $TMP/$RESULTS/nocoev
      cp -a $pair $TMP/$RESULTS/nocoev
      echo $SUMMARY >> $TMP/$RESULTS/coev_inter_all.tsv
    elif [ "$numPairs" -gt 0 ]; then
      echo -e "[\e[92mCOEVOL\e[39m] Copying pair with co-evolution: $pair"
      mkdir -p $TMP/$RESULTS/coev/$folder
      cp -a $pair $TMP/$RESULTS/coev/$folder
      echo $SUMMARY >> $TMP/$RESULTS/coev_inter.tsv
      echo $SUMMARY >> $TMP/$RESULTS/coev_inter_all.tsv
    else
      echo -e "Something went wrong for $pair ... Check!"
    fi
  done <<< $(echo "$SUMMARY")
}

# Cleanup results to make them easily parsable. Leave protein names
# and their UniProt identifiers in coev_inter.tsv. For *.out files,
# leave only the info of residues positions of the co-evolving pairs.
results_cleanup() {
  local coevPair="${1}"
  SUMMARY=$(sed -n '2p' $coevPair/coev_inter.csv)
  while read -r Seq1 Seq2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
    cp $coevPair/${Seq1}_${Seq2}.out $coevPair/$coevPair.clean
    echo "Clean up $coevPair.out"
    sed -i -n '/Coevolving Pairs of amino acid sites/,/Overlapping groups of coevolving residues/p' $coevPair/$coevPair.clean
    sed -i '1,5d' $coevPair/$coevPair.clean
    sed -i 's/Overlapping groups of coevolving residues//' $coevPair/$coevPair.clean
    sed -i '/^$/d' $coevPair/$coevPair.clean
    sed -i "s:\t\t:\t:g" $coevPair/$coevPair.clean
    sed -i "s/(/ /g" $coevPair/$coevPair.clean
    sed -i "s/)/ /g" $coevPair/$coevPair.clean
    sed -i "s/	/ /g" $coevPair/$coevPair.clean
    sed -i "s/  / /g" $coevPair/$coevPair.clean
    sed -i "s/^/$Seq1 $Seq2 /" $coevPair/$coevPair.clean
  done <<< $(echo "$SUMMARY")
}

columns_stats(){
  local coevPair="${1}"
  while read -r msa1 msa2 colA realA colB realB meanA meanB corr boot pvalA pvalB pMean corr1 corr2 ; do
    mkdir -p $coevPair/columnStats/$colA-$colB/
    extractalign \
     -sequence $coevPair/msa/$msa1 \
     -sbegin1 $colA \
     -send1 $colA \
     -sprotein1 \
     -osformat text \
     -osextension txt \
     -osname $coevPair/columnStats/$colA-$colB/$colA-$msa1 \
     -auto

    #grep ">" $coevPair/msa/$msa1 | sed "s:>::g" > $coevPair/columnStats/$colA-$colB/$colA-$msa1.species
    
    extractalign \
     -sequence $coevPair/msa/$msa2 \
     -sbegin1 $colB \
     -send1 $colB \
     -sprotein1 \
     -osformat text \
     -osextension txt \
     -osname $coevPair/columnStats/$colA-$colB/$colB-$msa2 \
     -auto

    #grep ">" $coevPair/msa/$msa2 | sed "s:>::g" > $coevPair/columnStats/$colA-$colB/$colB-$msa2.species

    #echo -e "Extracting MSA columns for $coevPair: $colA $colB"

    columnLenghtA=$(cat $coevPair/columnStats/$colA-$colB/$colA-$msa1.txt | wc -l)
    columnGapsA=$(grep "-" $coevPair/columnStats/$colA-$colB/$colA-$msa1.txt | wc -l)
    percentageA=$(echo "1 - ${columnGapsA}/${columnLenghtA}" | bc -l)
    roundPerGapA=$(printf "%1.5f" $percentageA)

    columnLenghtA=$(grep -v "-" $coevPair/columnStats/$colA-$colB/$colA-$msa1.txt | wc -l)
    columnUniqueA=$(grep -v "-" $coevPair/columnStats/$colA-$colB/$colA-$msa1.txt | sort | uniq -c | sort -n -r | sed -n 1p | awk '{ print $1 }')
    diversityResA=$(echo "1 - ${columnUniqueA}/${columnLenghtA}" | bc -l)
    roundDivResA=$(printf "%1.5f" $diversityResA)

    columnLenghtB=$(cat $coevPair/columnStats/$colA-$colB/$colB-$msa2.txt | wc -l)
    columnGapsB=$(grep "-" $coevPair/columnStats/$colA-$colB/$colB-$msa2.txt | wc -l)
    percentageB=$(echo "1 - ${columnGapsB}/${columnLenghtB}" | bc -l)
    roundPerGapB=$(printf "%1.5f" $percentageB)

    columnLenghtB=$(grep -v "-" $coevPair/columnStats/$colA-$colB/$colB-$msa2.txt | wc -l)
    columnUniqueB=$(grep -v "-" $coevPair/columnStats/$colA-$colB/$colB-$msa2.txt | sort | uniq -c | sort -n -r | sed -n 1p | awk '{ print $1 }')
    diversityResB=$(echo "1 - ${columnUniqueB}/${columnLenghtB}" | bc -l)
    roundDivResB=$(printf "%1.5f" $diversityResB)

    corrdec=$(printf "%1.10f" $corr)
    corr1dec=$(printf "%1.10f" $corr1)
    corr2dec=$(printf "%1.10f" $corr2)
    pvalAdec=$(printf "%1.10f" $pvalA)
    pvalBdec=$(printf "%1.10f" $pvalB)
    pMeandec=$(printf "%1.10f" $pMean)

    #echo -e "Calculating statistics for $coevPair: $colA $colB"
    
    echo "$msa1 $msa2 $colA $realA $colB $realB $meanA $meanB $corr1dec $corr2dec $corrdec $boot $roundPerGapA $roundPerGapB $roundDivResA $roundDivResB $pvalAdec $pvalBdec $pMeandec $coevPair" >> $coevPair/columnStatistics.tsv
    echo "$msa1-$msa2-$colA..$colB,$corr1dec,$corr2dec,$corrdec,$boot,$roundPerGapA,$roundPerGapB,$roundDivResA,$roundDivResB,$pvalAdec,$pvalBdec,$pMeandec" >> $TMP/$RESULTS/columnStatistics.csv

    # Remove pairs that do not pass the threshold of gaps, diversity, p-value and bootstrap.
    if (( $(echo "$roundPerGapA >= $RESGAPS" |bc -l) && \
          $(echo "$roundPerGapB >= $RESGAPS" |bc -l) && \
          $(echo "$roundDivResA >= $RESIDEN" |bc -l) && \
          $(echo "$roundDivResB >= $RESIDEN" |bc -l) && \
	  $(echo "$corr1dec > 0" |bc -l) && \
	  $(echo "$corr2dec > 0" |bc -l) && \
          $(echo "$boot >= $RESBOOT" |bc -l) && \
          $(echo "$pMeandec <= $PVALUE" |bc -l) )); then
      echo "$msa1 $msa2 $colA $realA $colB $realB $meanA $meanB $corr1dec $corr2dec $corrdec $boot $pvalAdec $pvalBdec $pMeandec" >> $coevPair/$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS
      echo -e "[\e[92mPASS\e[39m] Filter coev pair in $coevPair: $colA $colB"
    else
      echo "$msa1 $msa2 $colA $realA $colB $realB $meanA $meanB $corr1dec $corr2dec $corrdec $boot $pvalAdec $pvalBdec $pMeandec" >> $coevPair/$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-SKIP
      echo -e "[\e[37mSKIP\e[39m] Filter coev pair in $coevPair: $colA $colB"
    fi

  done < $coevPair/$coevPair.clean
}

mean_pval() {
  local coevPair="${1}"
  if find $coevPair/ -name "$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS" | read ; then
    BESTROW=$( sort --key 15 --numeric-sort $coevPair/$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS | head -n 1 )
    #BESTPVAL=$( awk 'BEGIN {bestp = 0} {if ($15<bestp) bestp=$15} END {print bestp}' $coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS)
    mpval=$( awk '{ meanp += $15 } END { print meanp/NR }' $coevPair/$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS )
    MEANPVAL=$(printf "%1.10f" $mpval)
    MEANBOOT=$( awk '{ meanb += $12 } END { print meanb/NR }' $coevPair/$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS )
    PAIRSNUM=$( wc -l $coevPair/$coevPair.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}-PASS | awk '{ print $1 }' )
    echo -e "[\e[95mCORR\e[39m] highest and mean values in $coevPair"
    echo "$BESTROW $MEANPVAL $MEANBOOT $PAIRSNUM" >> $TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}
    
  else
    echo -e "[\e[37mSKIP\e[39m] Good correlation not found for $coevPair"
  fi
}

summary_cleanup(){
  sed -i "s/_/ /g" $TMP/$RESULTS/coev_inter.tsv
  sed -i "s/\.fa/ /g" $TMP/$RESULTS/coev_inter.tsv
  sed -i "s/  / /g" $TMP/$RESULTS/coev_inter.tsv
  sed -i "s/_/ /g" $TMP/$RESULTS/coev_inter_all.tsv
  sed -i "s/\.fa/ /g" $TMP/$RESULTS/coev_inter_all.tsv
  sed -i "s/  / /g" $TMP/$RESULTS/coev_inter_all.tsv
  sed -i "s|_| |g" $TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}
  sed -i "s/\.fa/ /g" $TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}

  while read -r idxml namexml attribute ; do
    sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/coev_inter.tsv
    sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/coev_inter_all.tsv
    sed -i "s/$idxml/$namexml $idxml/g" $TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}
  done < $CWD/$PROTEIN

  sed -i "1i Protein1 UniProtID1 Protein2 UniProtID2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef" $TMP/$RESULTS/coev_inter.tsv
  sed -i "1i Protein1 UniProtID1 Protein2 UniProtID2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef" $TMP/$RESULTS/coev_inter_all.tsv
  sed -i "1i Pair,CorrelA,Correl2,Correl_mean,Bootstrap,gaps_A,gaps_B,identity_A,identity_B,PvalA,PvalB,Pval_mean" $TMP/$RESULTS/columnStatistics.csv
  sed -i "1i Protein1 UniProtID1 Protein2 UniProtID2 MSA1 real1 MSA2 real2 mean1 mean2 Correl1 Correl2 Correl Boot P-value1 P-value2 P-value P-all Boot-all Pairs-all" $TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}

  echo -e "\nPairs where CAPS run failed are placed in:\n\e[96m$TMP/$RESULTS/failed\e[39m"
  echo -e "\nPairs with no coevolving sites are placed in:\n\e[96m$TMP/$RESULTS/nocoev\e[39m"
  echo -e "\nPairs with coevolving sites are copied to:\n\e[96m$TMP/$RESULTS/coev\e[39m"
  echo -e "\nCleaned *.out files in:\n\e[96m$TMP/$RESULTS/coev\e[39m"
  echo -e "\nSummary of coev_inter.csv collected in:\n\e[96m$TMP/$RESULTS/coev_inter.tsv\e[39m"
}
