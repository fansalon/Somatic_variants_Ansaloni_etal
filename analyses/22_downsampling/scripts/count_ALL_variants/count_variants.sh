#! /bin/bash

set -e



############################################################################################
vcfdir=/projects/somatic_mouse_tissues/ONT/sniffles/02_variant_calling/output/*/*vcf*
out=/projects/somatic_mouse_tissues/ONT/downsampling/output/count_variants_allsamples/variant_count_ALL.txt
############################################################################################

echo "Counting variants in the vcf files in ${vcfdir}"




echo -e "sample\tVariants\tSomatic" > $out

for file in $vcfdir
do
  echo $file
  name=$(basename $file | awk -F ".vcf" '{ print $1 }')
  n_variant=$(zcat $file | grep -v "#" | wc -l)
  n_somatic=$(zcat $file | grep -v "#" | grep -w MOSAIC | wc -l)
  echo -e "${name}\t${n_variant}\t${n_somatic}" >> $out
done
