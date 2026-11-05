#! /bin/bash

set -e

###############################################################################
mergeR=scripts/merge_res.R
parsePy=scripts/parse_GT_tldr.py
refs=../../zenodo/input_data/teref.mouse.fa.fai
###############################################################################

echo "################################################################################"
echo "Starting:"


##### 1 - select only PASS, with TSD and not non-reference
# this outputs a file called tldr_res_filtered.txt
echo "Selecting PASS, with TSD and not non-reference MEI"
Rscript $mergeR
nro1=$(wc -l input/tldr_res.txt | awk '{ print $1 }')
nro2=$(wc -l tldr_res_filtered.txt | awk '{ print $1 }')
echo -e "\t${nro2}/${nro1} MEI selected"


##### 2 - convert file to bed, parsing GT, REF and ALT, defining somatic events
echo "Selecting somatic MEIs"
python3 $parsePy --input tldr_res_filtered.txt --refs $refs
nro3=$(cat tldr_res_filtered.txt.bed | cut -f 16 | grep -c True)
echo -e "\t${nro3}/${nro2} somatic MEI selected"



echo Done
echo "################################################################################"
