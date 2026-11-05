#########################################################################################################################################################
### 12.
### Overlap tissue-specific somatic private BNDs with TEs
#########################################################################################################################################################
rm(list = ls())
t00=Sys.time()
#########################################################################################################################################################
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(ChIPseeker))
suppressMessages(library(BiocParallel))
suppressMessages(library(parallel))
suppressMessages(library(rlist))
suppressMessages(library(regioneR))
#########################################################################################################################################################
source("../functions/functions_nanopore.R")
#########################################################################################################################################################
path_to_rmsk='../../zenodo/input_data/rmsk.mm39.txt.filtered.bed'
bedtools="/opt/homebrew/bin/bedtools"
txdb_path='../../zenodo/input_data/txdb.db'
#########################################################################################################################################################
bootstrap=1000
seed_f=110322
#########################################################################################################################################################
# txdb
txdb <- loadDb(txdb_path)
# Somatic anno
obj <- readRDS("../02_get_somatic/somatic.Rds")
comparison <- obj$tissue_specificity
comparison$Type <- gsub("_vs_","-vs-",comparison$Type)
# Somatic SVs
df <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
#########################################################################################################################################################

cat("\n##################################################################################\n")
cat(paste0("Starting analysis:\n\tBootstrap: ",bootstrap,"\n"))
cat("Selecting SVs of interest:\n")

# get Somatic_tissue_private
df1 <- subset(df,somatic.anno=="Somatic_tissue_private")
cat(paste0("\tSelected ", nrow(df1)," somatic tissue-specific private variants\n"))
# get BND
df2 <- subset(df1,SVTYPE=="BND")
cat(paste0("\tSelected ", nrow(df2)," somatic tissue-specific private BNDs\n"))

# add tissue comparison column
tissue_specific2_expand <- addComparisonType(df2,comparison)
table(tissue_specific2_expand$Type)
cat("\n")
#########################################################################################################################################################


cat("Converting to bed\n")

# convert to bed
bed1 <- tissue_specific2_expand[,c("chrom","start","end","ID","Type")]
bed1$type <- "bp1"
bed2 <- tissue_specific2_expand[,c("breakpoint2","ID","Type")]
bed2$chrom <- str_split_fixed(bed2$breakpoint2,":",2)[,1]
bed2$start <- as.numeric(str_split_fixed(bed2$breakpoint2,":",2)[,2])
bed2$end <- bed2$start+1
bed2 <- bed2[,c("chrom","start","end","ID","Type")]
bed2$type <- "bp2"
bed <- rbind(bed1,bed2)



cat("Overlapping with TEs\n")

# overlap TEs - REAL
te_overlap <- wrap_overlapBND_TEs(bed,path_to_rmsk,bedtools,i=1,TEsplit=2)



cat("Shuffling analysis\n")

# overlap TEs - SHUFFLE
t0=Sys.time()
te_overlap_shuf <- wrap_overlapBND_TEs_shuf(bed,txdb,path_to_rmsk,bedtools,TEsplit=2,bootstrap)
t1=Sys.time()
print(t1-t0)

# bootstrap    2: 1'
# bootstrap   10: 3.2'
# bootstrap  100: 30'
# bootstrap 1000:
#########################################################################################################################################################


cat("Calculating O/E statistics aggregating <1% TE pairs\n")


# merge obs and exp
te_overlap$ID <-paste0(te_overlap$TE,"@",te_overlap$sample)
te_overlap$TE <- NULL
te_overlap$sample <- NULL
te_overlap_shuf$ID <-paste0(te_overlap_shuf$TE,"@",te_overlap_shuf$sample)
te_overlap_shuf$TE <- NULL
te_overlap_shuf$sample <- NULL
merged <- merge(te_overlap,te_overlap_shuf,by="ID",all.x=T)
merged[is.na(merged)] = 0
names(merged)[names(merged)=="Freq"] <- "Obs"

# select TEs with both Obs and Exp >0
merged <- subset(merged,Obs > 0 & Expected > 0)


# split comparison information and TE annotation
merged$comparison <- str_split_fixed(merged$ID,"@",2)[,2]
merged$TE <- str_split_fixed(merged$ID,"@",2)[,1]
merged$ID <- NULL

# add totals
tots <- aggregate(merged$Obs,by=list(merged$comparison),sum)
names(tots) <- c("comparison","Tot")
merged <- merge(merged,tots,by="comparison",all.x=T)

# calc %
merged$Obs_perc <- merged$Obs / merged$Tot * 100
merged$Exp_perc <- merged$Expected / merged$Tot * 100

# classify <1% as 'other'
merged.orig <- merged
merged$TE <- ifelse(merged$Obs_perc<1,"Other",as.character(merged$TE))

# aggregate
merged$var <- merged$sd^2
merged_a <- aggregate(merged$Obs,by=list(merged$TE,merged$comparison),sum)
names(merged_a) <- c("TE","comparison","Obs")
merged_a$Obs_perc <- aggregate(merged$Obs_perc,by=list(merged$TE,merged$comparison),sum)[,"x"]
merged_a$Exp <- aggregate(merged$Expected,by=list(merged$TE,merged$comparison),sum)[,"x"]
merged_a$Exp_perc <- aggregate(merged$Exp_perc,by=list(merged$TE,merged$comparison),sum)[,"x"]
merged_a$var <- aggregate(merged$var,by=list(merged$TE,merged$comparison),sum)[,"x"]
merged_a$sd <- sqrt(merged_a$var) 
merged_a$var <- NULL
merged_a$bootstrap <- bootstrap

# add totals
merged_a <- merge(merged_a,tots,by="comparison",all.x=T)
merged_a$sd_perc <- merged_a$sd / merged_a$Tot * 100

# calc zscore, pval and fdr
merged_a2 <- zscoreCalc(merged_a)
#########################################################################################################################################################


cat("Calculating O/E statistics on all the TE pairs\n")


# calc zscore, pval and fdr - non aggregated file
colnames(merged.orig)[colnames(merged.orig)=="Expected"] <- "Exp"
merged.orig$sd_perc <- merged.orig$sd / merged.orig$Tot * 100
merged.orig2 <- zscoreCalc(merged.orig)
merged.orig2 <- merged.orig2[,colnames(merged_a2)]


#########################################################################################################################################################

cat("Writing to out\n")

for (j in 1:ncol(merged.orig2)) {
  if (class(merged.orig2[,j])=="numeric") {
    if (colnames(merged.orig2)[j] != "fdr" & colnames(merged.orig2)[j] != "pvalue") {
      merged.orig2[,j] <- round(merged.orig2[,j],2)
      merged_a2[,j] <- round(merged_a2[,j],2)
    } else {
      merged.orig2[,j] <- formatC(merged.orig2[,j],digits = 2,format = 'e')
      merged_a2[,j] <- formatC(merged_a2[,j],digits = 2,format = 'e')
    }
  }
}



write.table(merged_a2,paste0("BND_overlapTEs_aggr_res_boot",bootstrap,".txt"),row.names = F,sep='\t',quote=F)
write.table(merged.orig2,paste0("BND_overlapTEs_res_boot",bootstrap,".txt"),row.names = F,sep='\t',quote=F)

t11=Sys.time()
print(t11-t00)

cat("Done\n##################################################################################\n")















