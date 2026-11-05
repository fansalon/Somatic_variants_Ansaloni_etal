#! /bin/bash

set -e


if [ $# -ne 5 ]
then
  echo "missing args"
  echo "usage: <bam> <fasta genome> <cpu> <outdir> <tandem repeats from findTandemRepeats.py>"
  exit
fi


### Input param
bam=$1
genome=$2
cpu=$3
outdir=$4
tr=$5


# create names and outdir name
name=$(basename $bam | awk -F ".bam" '{ print $1 }')
outdir=${outdir}/${name}
if [ ! -d $outdir ]; then mkdir $outdir; fi
cd $outdir


# define vcf and snf file name
vcf=${name}.vcf
snf=${name}.snf



###
# 1. Launch sniffles
date
echo "Launching sniffles on $bam..."

singularity exec --bind /projects/somatic_mouse_tissues/ONT/ /projects/somatic_mouse_tissues/ONT/sniffles/src/sniffles2.2.simg sniffles --input $bam --reference $genome --threads $cpu --output-rnames --minsupport 3 --minsvlen 20 --sample-id $name --snf $snf --vcf $vcf --mosaic --mosaic-af-min 0.03 --mosaic-include-germline --mosaic-qc-invdup-min-length 250 --tandem-repeats $tr


echo done
date


