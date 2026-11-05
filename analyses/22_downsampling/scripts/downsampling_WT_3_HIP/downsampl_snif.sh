#!/bin/bash

if [ $# -ne 3 ]
then
  echo "ERROR"
  echo "args: <bam> <number of the boostrap> <output dir>"
  exit
fi

cd $PBS_O_WORKDIR

################################################################
module load samtools
genome=/projects/somatic_mouse_tissues/ONT/sniffles/genome/GRCm39.primary_assembly.genome.fa
tr=/projects/somatic_mouse_tissues/ONT/sniffles/01_get_tandem_rep/tandem_rep_GRCm39.bed
sniffles_exe="singularity exec --bind /projects/somatic_mouse_tissues/ONT/ /projects/somatic_mouse_tissues/ONT/sniffles/src/sniffles2.2.simg sniffles"
################################################################


################################################################
bam=$1
bootstrap_n=$2
outdir=$3
out=${outdir}/${bootstrap_n}
mkdir $out
cd $out
echo -e "Downsample\tVariants\tSomatic\tMapped_reads" > variant_count.txt
################################################################



for i in $(seq 1 9)
do 
  date
  echo -e "\nDownsample: 0.${i}"

  # generate dir
  mkdir down_0.${i}
  cd down_0.${i}

  # samtools
  echo -e "\tgenerating bam.."
  mkdir bam
  cd bam
  samtools view -@ 20 --subsample "0.${i}" --subsample-seed 110322${bootstrap_n} -o subam_.${i}.bam $bam
  #samtools view -@ 20 -subsample 110322${bootstrap_n}".${i}" -o subam_.${i}.bam $bam
  samtools index -@ 20 subam_.${i}.bam
  samtools flagstat -@ 20 subam_.${i}.bam > subam_.${i}.flagstat
  mapped=$(grep "+ 0 mapped" subam_.${i}.flagstat | awk '{ print $1 }')
  newbam=$PWD/subam_.${i}.bam
  cd ../

  # call peaks
  date
  echo -e "\tcounting variants.."
  mkdir variants
  cd variants
  $sniffles_exe --input $newbam --reference $genome --threads 20 --minsupport 3 --minsvlen 20 --vcf variants_.${i}.vcf \
	  --mosaic --mosaic-af-min 0.03 --mosaic-include-germline --mosaic-qc-invdup-min-length 250 --tandem-repeats $tr --quiet 

  # count variants
  n_variant=$(cat variants_.${i}.vcf | grep -v "#" | wc -l)
  n_somatic=$(cat variants_.${i}.vcf| grep -v "#" | grep -w MOSAIC | wc -l)

  # gzip vcf
  gzip variants_.${i}.vcf

  # write to out
  cd ../
  echo -e ".${i}\t${n_variant}\t${n_somatic}\t${mapped}" >> ../variant_count.txt

  # remove bam
  date
  echo -e "\tremoving bams.."
  rm bam/*bam bam/*bai

  cd ../
done



################################################################
# full bam
mkdir down_1
cd down_1
date
echo -e "\nFull bam"
mkdir bam
cd bam

# flagstat
echo -e "\tgenerating bam.. (skipped)"
samtools flagstat -@ 20 ${bam} > subam_1.flagstat
mapped=$(grep "+ 0 mapped" subam_1.flagstat | awk '{ print $1 }')
cd ../

# variants
mkdir variants
cd variants
echo -e "\tcounting variants.."
$sniffles_exe --input $bam --reference $genome --threads 20 --minsupport 3 --minsvlen 20 --vcf variants_1.vcf \
	--mosaic --mosaic-af-min 0.03 --mosaic-include-germline --mosaic-qc-invdup-min-length 250 --tandem-repeats $tr --quiet 

# count variants
n_variant=$(cat variants_1.vcf | grep -v "#" | wc -l)
n_somatic=$(cat variants_1.vcf| grep -v "#" | grep -w MOSAIC | wc -l)
gzip variants_1.vcf
cd ..
echo -e "1\t${n_variant}\t${n_somatic}\t${mapped}" >> ../variant_count.txt
################################################################


date
echo "ALL Done"

