#!/bin/bash

# Generate nodes list
node_gen_list() {
  sed 1d $TMP/$RESULTS/coev_inter.tsv | while read -r Protein1 UniProtID1 Protein2 UniProtID2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
    echo $UniProtID1 $Protein1 >> $TMP/$RESULTS/nodes.txt
    echo $UniProtID2 $Protein2 >> $TMP/$RESULTS/nodes.txt
  done

  sort $TMP/$RESULTS/nodes.txt | uniq > $TMP/$RESULTS/nodesUnique.txt
  echo -e "Generated nodes list."
}

make_coev_inter_xml() {
# Prepare the XML header
cat << EOF >> $TMP/$RESULTS/CAPS.coev_inter.xml
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<graph label="1478" directed="0" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">
  <att name="networkData">
    <rdf:RDF>
      <rdf:Description rdf:about="http://www.cytoscape.org/">
        <dc:type>Protein-Protein Interaction</dc:type>
        <dc:description>N/A</dc:description>
        <dc:identifier>N/A</dc:identifier>
        <dc:date>$DATE</dc:date>
        <dc:title>CAPS</dc:title>
        <dc:source>http://www.cytoscape.org/</dc:source>
        <dc:format>Cytoscape-XGMML</dc:format>
      </rdf:Description>
    </rdf:RDF>
  </att>
  <att name="name" value="coev_inter: Alpha: ${ALPHA}; MSA: ${MSAMETHOD}; Tree: ${TREESCAPS}" type="string" cy:type="String"/>
  <att name="Publication" value="Parallelized CAPS2." type="string" cy:type="String"/>
  <att name="Dataset Name" value="APEX biotinylation" type="string" cy:type="String"/>
  <att name="Dataset URL" value="https://www.github.com/mattilalab" type="string" cy:type="String"/>
EOF

# Add dummy nodes and edges
cat << EOF >> $TMP/$RESULTS/CAPS.coev_inter.xml
  <node id="0000XXXX" label="DummyNode1">
    <att name="shared name" value="DummyNode1" type="string" cy:type="String"/>
    <att name="name" value="DummyNode1" type="string" cy:type="String"/>
  </node>
  <node id="1111ZZZZ" label="DummyNode2">
    <att name="shared name" value="DummyNode2" type="string" cy:type="String"/>
    <att name="name" value="DummyNode2" type="string" cy:type="String"/>
  </node>
  <edge id="0000XXXX-pos-1111ZZZZ" label="DummyNode1 (pos_coev) DummyNode2" source="0000XXXX" target="1111ZZZZ" cy:directed="0">
    <att name="shared name" value="DummyNode1 (pos_coev) DummyNode2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="DummyNode1 (pos_coev) DummyNode2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="numPairs" value="1" type="real" cy:type="Double"/>
    <att name="thresholdR" value="0" type="real" cy:type="Double"/>
    <att name="averageR" value="0" type="real" cy:type="Double"/>
    <att name="averageSigR" value="1" type="real" cy:type="Double"/>
  </edge>
  <edge id="0000XXXX-neg-1111ZZZZ" label="DummyNode1 (neg_coev) DummyNode2" source="0000XXXX" target="1111ZZZZ" cy:directed="0">
    <att name="shared name" value="DummyNode1 (neg_coev) DummyNode2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="DummyNode1 (neg_coev) DummyNode2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="numPairs" value="1" type="real" cy:type="Double"/>
    <att name="thresholdR" value="0" type="real" cy:type="Double"/>
    <att name="averageR" value="0" type="real" cy:type="Double"/>
    <att name="averageSigR" value="-1" type="real" cy:type="Double"/>
  </edge>
EOF

# Add nodes.
while read -r idxml namexml ; do
echo -e "Adding node $idxml $namexml in CAPS.coev_inter.xml"
cat << EOF >> $TMP/$RESULTS/CAPS.coev_inter.xml
  <node id="$idxml" label="$namexml">
    <att name="shared name" value="$namexml" type="string" cy:type="String"/>
    <att name="name" value="$namexml" type="string" cy:type="String"/>
  </node>
EOF
done < $TMP/$RESULTS/nodesUnique.txt

# Add edges
sed 1d $TMP/$RESULTS/coev_inter.tsv | while read -r Protein1 UniProtID1 Protein2 UniProtID2 numPairs totalComp CutOff thresholdR averageR averageSigR tree1length tree2length gapThreshold bootCutOff DistanceCoef ; do
echo -e "Adding edge: $UniProtID1-$UniProtID2 in CAPS.coev_inter.xml"
cat << EOF >> $TMP/$RESULTS/CAPS.coev_inter.xml
  <edge id="$UniProtID1-$UniProtID2" label="$Protein1 (coev) $Protein2" source="$UniProtID1" target="$UniProtID2" cy:directed="0">
    <att name="shared name" value="$Protein1 (coev) $Protein2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="$Protein1 (coev) $Protein2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="numPairs" value="$numPairs" type="real" cy:type="Double"/>
    <att name="thresholdR" value="$thresholdR" type="real" cy:type="Double"/>
    <att name="averageR" value="$averageR" type="real" cy:type="Double"/>
    <att name="averageSigR" value="$averageSigR" type="real" cy:type="Double"/>
  </edge>
EOF
done
echo "</graph>" >> $TMP/$RESULTS/CAPS.coev_inter.xml
}

make_pos_corr_xml() {
if [ -f "$TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}" ]; then
# Prepare the XML header
cat << EOF >> $TMP/$RESULTS/CAPS.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}.xml
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<graph label="1478" directed="0" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">
  <att name="networkData">
    <rdf:RDF>
      <rdf:Description rdf:about="http://www.cytoscape.org/">
        <dc:type>Protein-Protein Interaction</dc:type>
        <dc:description>N/A</dc:description>
        <dc:identifier>N/A</dc:identifier>
        <dc:date>$DATE</dc:date>
        <dc:title>CAPS</dc:title>
        <dc:source>http://www.cytoscape.org/</dc:source>
        <dc:format>Cytoscape-XGMML</dc:format>
      </rdf:Description>
    </rdf:RDF>
  </att>
  <att name="name" value="PosCorr; MSA: ${MSAMETHOD}; Tree: ${TREESCAPS}" type="string" cy:type="String"/>
  <att name="Publication" value="Parallelized CAPS2. Petrov P, Sustar V, Awoniyi L, Toivonen M, Mattila P" type="string" cy:type="String"/>
  <att name="Dataset Name" value="APEX biotinylation" type="string" cy:type="String"/>
  <att name="Dataset URL" value="https://www.github.com/mattilalab/Petrov-Sustar_et_al" type="string" cy:type="String"/>
EOF

# Add nodes.
sed 1d $TMP/tsv/proteinsFound.tsv | while read -r idxml namexml attribute; do
echo -e "Adding node $idxml $namexml $attribute in $TMP/$RESULTS/CAPS.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}.xml"
cat << EOF >> $TMP/$RESULTS/CAPS.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}.xml
  <node id="$idxml" label="$namexml">
    <att name="shared name" value="$namexml" type="string" cy:type="String"/>
    <att name="name" value="$namexml" type="string" cy:type="String"/>
    <att name="attribute" value="$attribute" type="string" cy:type="String"/>
  </node>
EOF
done

# Add edges
sed 1d $TMP/$RESULTS/stats.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE} | while read -r Protein1 UniProtID1 Protein2 UniProtID2 MSA1 real1 MSA2 real2 mean1 mean2 Correl1 Correl2 Correl Boot Pvalue1 Pvalue2 Pvalue Pall Bootall Pairsall ; do
echo -e "Adding edge: $UniProtID1-$UniProtID2 in $TMP/$RESULTS/CAPS.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}.xml"
cat << EOF >> $TMP/$RESULTS/CAPS.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}.xml
  <edge id="$UniProtID1-$UniProtID2" label="$Protein1 (pos_coev) $Protein2" source="$UniProtID1" target="$UniProtID2" cy:directed="0">
    <att name="shared name" value="$Protein1 (pos_coev) $Protein2" type="string" cy:type="String"/>
    <att name="shared interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="name" value="$Protein1 (pos_coev) $Protein2" type="string" cy:type="String"/>
    <att name="interaction" value="Coevolution" type="string" cy:type="String"/>
    <att name="Best P-value" value="$Pvalue" type="real" cy:type="Double"/>
    <att name="Best pair bootstrap" value="$Boot" type="real" cy:type="Double"/>
    <att name="Mean P-value" value="$Pall" type="real" cy:type="Double"/>
    <att name="Mean bootstrap" value="$Bootall" type="real" cy:type="Double"/>
    <att name="Number of Pairs" value="$Pairsall" type="real" cy:type="Double"/>
  </edge>
EOF
done
echo "</graph>" >> $TMP/$RESULTS/CAPS.G${RESGAPS}-I${RESIDEN}-B${RESBOOT}-P${PVALUE}.xml
else
  echo -e "No codependence found"
fi
}
