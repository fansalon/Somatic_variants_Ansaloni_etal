###########################################################################################################
### 04.
### This script gets all the SV and annotated them compared to genomic features also 
### performing bootstrapping.
###########################################################################################################
rm(list = ls())
seed_f <- "11322"
set.seed(seed_f)
###########################################################################################################
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(BiocParallel))
suppressMessages(library(rlist))
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(regioneR))
suppressMessages(library(stringr))
suppressMessages(library(parallel))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################
priorities <- c("Promoter","utr5","utr3","Exon","Intron","Downstream")
bootstrap <- 1000
###########################################################################################################

cat("\n##################################################################################\n")
cat("Loading:\n\t../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed\n")

# Load vcf
vcfbed <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
# select only variants of interest
cat("\tselecting:\n\t\tSomatic_common_private\n\t\tSomatic_common_recurrent\n\t\tSomatic_tissue_private\n")
vcfbed <- subset(vcfbed,somatic.anno=="Somatic_common_private"|
                   somatic.anno=="Somatic_common_recurrent"|
                   somatic.anno=="Somatic_tissue_private")


# Load txdb
cat("Loading:\n\ttxdb.db\n")
txdb <- loadDb("../../zenodo/input_data/txdb.db")

# Load somatic obj
cat("Loading:\n\t../02_get_somatic/somatic.Rds\n")
obj <- readRDS("../02_get_somatic/somatic.Rds")

# get genomic features coords in txdb
cat("Extracting features from txdb\n")
feature_Grl <- extractFeatureFromTxDb(TxDb=txdb,upstreamTSS=3000,downstreamTSS=3000,downstreamTTS=300)
###########################################################################################################


### Annotate: common recurrent and private DEL and INS and tissue specific private BND
### For convenience this is done in a wrapper function, as the very same operation will
### be repeated for variants in simple rep and normal regions
cat("#################################################################################\n")
cat("Annotating ALL the variants\n")
ff <- wrap_GenomicFeat(vcfbed,sv_common=c("DEL","INS"),sv_tissue="BND")
write.table(ff,paste0("final_anno_table_ALL_",bootstrap,".txt"),row.names = F,sep='\t',quote=F)



### The same for variants in simple repeats
cat("#################################################################################\n")
cat("Annotating only the variants in simple repeat regions\n")
vcfbed_SR <- subset(vcfbed,final.anno=="simple_rep")
ff_SR <- wrap_GenomicFeat(vcfbed_SR,sv_common=c("DEL","INS"),sv_tissue="BND")
write.table(ff_SR,paste0("final_anno_table_SimpleRep_",bootstrap,".txt"),row.names = F,sep='\t',quote=F)


### The same for variants NON in simple repeats
cat("#################################################################################\n")
cat("Annotating only the variants not in simple repeat regions\n")
vcfbed_nSR <- subset(vcfbed,final.anno!="simple_rep")
ff_nSR <- wrap_GenomicFeat(vcfbed_nSR,sv_common=c("DEL","INS"),sv_tissue="BND")
write.table(ff_nSR,paste0("final_anno_table_normal_",bootstrap,".txt"),row.names = F,sep='\t',quote=F)

cat("Done.\n")



