###########################################################################################################
### 14.
### This script annotates respect to genomic features the somatic MEI
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
gene_anno <- fread("../../zenodo/input_data/gencode.GRCm39.vM28.annotation.gtf.transcript.bed",
                   sep='\t',data.table = F)[,c(5,7)]
names(gene_anno) <- c("tx_name","gene_name")
###########################################################################################################


cat("\n##################################################################################\n")
cat("Loading:\n\t../13_tldr/tldr_res_filtered.txt.bed\n")

# Load vcf
vcfbed <- fread("../13_tldr/tldr_res_filtered.txt.bed",header=T,sep='\t',data.table = F)
# select only somatic variants
vcfbed <- subset(vcfbed,is.somatic==TRUE)
cat(paste0("\tselecting:\n\t\tsomatic MEIs (",nrow(vcfbed),")\n"))

# Load txdb
cat("Loading:\n\ttxdb.db\n")
txdb <- loadDb("../../zenodo/input_data/txdb.db")

# get genomic features coords in txdb
cat("Extracting features from txdb\n")
feature_Grl <- extractFeatureFromTxDb(TxDb=txdb,upstreamTSS=3000,downstreamTSS=3000,downstreamTTS=300)
###########################################################################################################


### Annotate
cat("#################################################################################\n")
cat("Annotating somatic MEIs\n")



############# Real
# select only INS breakpoint
sv_of_i <- vcfbed[,c("Chrom","Start","ID")]
sv_of_i$end <- sv_of_i$Start+1
sv_of_i <- sv_of_i[,c("Chrom","Start","end","ID")]
names(sv_of_i) <- c("chrom","start","end","ID")
sv_of_i$sample_comparison <- "sample01"

sv_of_i.ann <- annotateGenomicFeatures2(df=sv_of_i,
                                        priorities = priorities,
                                        feature_Grl,
                                        quiet = T,
                                        is.shuffle = F,
                                        txdb_raw=txdb)
sv_of_i.ann$id <- NULL
sv_of_i.ann$chrom <- NULL
sv_of_i.ann$start <- NULL
sv_of_i.ann$end <- NULL
sv_of_i.ann$sample_comparison <- NULL

# add gene annotation
sv_of_i.ann <- merge(sv_of_i.ann,gene_anno,by="tx_name",all.x=T,sort=F)




############# Shuffle
cat("Annotating shuffle somatic MEIs\n")
t0 <- Sys.time()
exp <- wrapShuffleCount2(df = sv_of_i,
                         txdb = txdb,
                         prior = priorities,
                         txdb_GR = feature_Grl,
                         bootstrap = bootstrap,
                         threads=5,
                         sample_id="sample_comparison",
                         is.shuffle = T,
                         txdb_raw = txdb)
t1 <- Sys.time()
t1-t0

# bootstrap 2    = 2.5''
# bootstrap 10   = 4''
# bootstrap 100  = 31''
# bootstrap 1000 = 3.4'


### calc zscore
obs <- as.data.frame(table(sv_of_i.ann$annotation))
names(obs) <- c("annotation","Obs")

# merge obs exp
merged <- merge(obs,exp,by='annotation',all.y=T)[,c("annotation","Obs","Exp","sd","bootstrap")]
merged[is.na(merged)] = 0
merged$Tot <- sum(merged$Obs)

# calc perc
merged$Obs_perc <- merged$Obs / merged$Tot * 100
merged$Exp_perc <- merged$Exp / merged$Tot * 100

# calc zscore, pvalue, fdr
merged <- zscoreCalc(merged)
# format
format1 <- c("Exp","sd","Obs_perc","Exp_perc","Enrichment","zscore")
for (col in format1) {
  merged[,col] <- round(merged[,col],2)
}
format2 <- c("pvalue","fdr")
for (col in format2) {
  merged[,col] <- formatC(merged[,col],digits = 2,format = 'e')
}



############# write to out
write.table(sv_of_i.ann,"genomic_feature_anno.txt",row.names = F,sep='\t',quote=F)
write.table(merged,"res_observed_expected.txt",row.names = F,sep='\t',quote=F)



cat("Done.\n")
cat("#################################################################################\n")


