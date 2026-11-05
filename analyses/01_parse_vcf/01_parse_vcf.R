###########################################################################################################
### 01.
### This script converts the original vcf file to a vcf.bed file
###########################################################################################################
rm(list = ls())
suppressMessages(library(stringr))
suppressMessages(library(BiocParallel))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(txdbmaker))
suppressMessages(library(regioneR))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################
### Please define your own
bedtools <- "/opt/homebrew/bin/bedtools"
###########################################################################################################
### Use at it is
inpvcf <- "../../zenodo/input_data/Table_S1.vcf"
gtf <- "../../zenodo/input_data/gencode.GRCm39.vM28.annotation.gtf"
genomeN <- "../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr"
gene_annoN <- "../../zenodo/input_data/gencode.GRCm39.vM28.annotation.gtf.transcript.bed"
rmsk <- "../../zenodo/input_data/rmsk.mm39.txt.filtered.bed"
simple_rep1 <- "../../zenodo/input_data/simpleRepeat.txt.bed"
simple_rep2 <- "../../zenodo/input_data/simpleRepeat.rmsk.bed"
priorities <- c("Promoter","utr5","utr3","Exon","Intron","Downstream")
###########################################################################################################


# Load vcf
cat("\n#############################################################\n")
cat("Loading vcf (Table_S1)\n")
vcf <- read.table(inpvcf,comment.char = "#",sep='\t',header = F)
colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                   "WT_2_CTX","WT_2_HIP","WT_2_KID","WT_3_CTX","WT_3_HIP","WT_3_KID",
                   "WT_4_CTX","WT_4_HIP","WT_4_KID","WT_5_CTX","WT_5_HIP","WT_5_KID","WT_6_CTX","WT_6_HIP","WT_6_KID")

# Remove unused columns
vcf$REF <- NULL
vcf$ID <- NULL
vcf$FORMAT <- NULL


### Get info from INFO col (precision, svtype, svlen, etc)
cat("Parsing INFO vcf field\n")
info_list <- str_split_fixed(vcf$INFO,";",12)
vcf$INFO <- NULL
vcf$PRECISION <- info_list[,1]
vcf$SVTYPE <- gsub("SVTYPE=","",info_list[,2])

# the format of the following is different for BND and all the other SVs
vcf$SVLEN <- as.numeric(ifelse(vcf$SVTYPE!="BND",gsub("SVLEN=","",info_list[,3]),NA))
vcf$SVLEN <- abs(vcf$SVLEN)
vcf$END <- as.numeric(ifelse(vcf$SVTYPE!="BND",gsub("END=","",info_list[,4]),vcf$POS+1))
vcf$SUPPORT <- as.numeric(ifelse(vcf$SVTYPE!="BND",gsub("SUPPORT=","",info_list[,5]),gsub("SUPPORT=","",info_list[,3])))
vcf$COVERAGE <- ifelse(vcf$SVTYPE!="BND",gsub("COVERAGE=","",info_list[,7]),gsub("COVERAGE=","",info_list[,5]))
vcf$STRAND <- ifelse(vcf$SVTYPE!="BND",gsub("STRAND=","",info_list[,8]),gsub("STRAND=","",info_list[,6]))
vcf$ALLELE_COUNT <- as.numeric(ifelse(vcf$SVTYPE!="BND",gsub("AC=","",info_list[,9]),gsub("AC=","",info_list[,7])))
vcf$STDEV_LEN <- as.numeric(ifelse(vcf$SVTYPE!="BND",gsub("STDEV_LEN=","",info_list[,10]),gsub("STDEV_LEN=","",info_list[,8])))
vcf$STDEV_POS <- as.numeric(ifelse(vcf$SVTYPE!="BND",gsub("STDEV_POS=","",info_list[,11]),gsub("STDEV_POS=","",info_list[,9])))
vcf$SUPP_VEC <- ifelse(vcf$SVTYPE!="BND",gsub("SUPP_VEC=","",info_list[,12]),gsub("SUPP_VEC=","",info_list[,10]))


### Get second breakpoint for BND
vcf$breakpoint2 <- ifelse(vcf$SVTYPE!="BND",NA,vcf$ALT)
# remove "[" or "]" from breakpoint
vcf$breakpoint2 <- ifelse(vcf$SVTYPE!="BND",NA,
                          ifelse(grepl("\\[",vcf$breakpoint2),str_split_fixed(vcf$ALT,"\\[",3)[,2],
                                 str_split_fixed(vcf$ALT,"\\]",3)[,2]))

# do a check
if (! table(is.na(subset(vcf,SVTYPE=="BND")[,"breakpoint2"]))["FALSE"]==nrow(subset(vcf,SVTYPE=="BND"))) {
  stop("ERROR: there are some BND not assigned to any breakpoint2")
}
vcf$ALT <- NULL


### Reorder samples --> all CTX, all HIP and all KID
ctx <- vcf[,grepl("CTX",colnames(vcf))]
hip <- vcf[,grepl("HIP",colnames(vcf))]
kid <- vcf[,grepl("KID",colnames(vcf))]
# remove columns
vcf <- vcf[,!grepl("CTX",colnames(vcf))]
vcf <- vcf[,!grepl("HIP",colnames(vcf))]
vcf <- vcf[,!grepl("KID",colnames(vcf))]
# add in correct order
vcf <- cbind(vcf,ctx,hip,kid)


### parse the info of the SV called in each single sample
cat("Parsing genotype and supporting read fields for each sample\n")
counter <- 0
for (col in colnames(vcf)) {
  if (grepl("WT",col)) {
    counter <- counter + 1
    #print(col)
    
    # get all the info, sample per sample (e.g. 0/1:60:36:35:Sniffles2.DEL.204AS0)
    info <- vcf[,col]
    
    # split them into dataframes
    info_df <- str_split_fixed(info,":",5)
    
    # get info of interest (genotype, genotype qual, no. reference reads, no. variant reads)
    info_df_GT <- data.frame(info_df[,1]) # genotype
    colnames(info_df_GT) <- paste0(col,"_GT")
    info_df_GQ <- data.frame(as.numeric(info_df[,2])) # genotype qual
    colnames(info_df_GQ) <- paste0(col,"_GQ")
    info_df_DR <- data.frame(as.numeric(info_df[,3])) # no. reference reads
    colnames(info_df_DR) <- paste0(col,"_REF")
    info_df_DV <- data.frame(as.numeric(info_df[,4])) # no. variant reads
    colnames(info_df_DV) <- paste0(col,"_ALT")
    
    # combine together the same type of info for the different samples
    if (counter==1) {
      df_GT <- info_df_GT # genotype
      df_GQ <- info_df_GQ # genotype qual
      df_DR <- info_df_DR # no. reference reads
      df_DV <- info_df_DV # no. variant reads
    } else {
      df_GT <- cbind(df_GT,info_df_GT) # genotype
      df_GQ <- cbind(df_GQ,info_df_GQ) # genotype qual
      df_DR <- cbind(df_DR,info_df_DR) # no. reference reads
      df_DV <- cbind(df_DV,info_df_DV) # no. variant reads
    }
    vcf[,col] <- NULL
  }
}


### combine together the values for the 18 samples (for genotype qual and reference reads)
vcf$GENOTYPE.QUAL_VEC <- apply( df_GQ[ ,  ] , 1 , paste , collapse = "|" )
vcf$REF_READS <- apply( df_DR[ ,  ] , 1 , paste , collapse = "|" )

### add columns referring to genotype of the SV and ALT supporting reads in each sample 
vcf <- cbind(vcf,df_GT,df_DV)


### Count number of samples with GT 0/0 in the 6 samples of each tissue
cat("Calculating metrics\n")
ctx <- c("WT_2_CTX_GT","WT_3_CTX_GT","WT_4_CTX_GT","WT_5_CTX_GT","WT_6_CTX_GT")
hip <- c("WT_2_HIP_GT","WT_3_HIP_GT","WT_4_HIP_GT","WT_5_HIP_GT","WT_6_HIP_GT")
kid <- c("WT_2_KID_GT","WT_3_KID_GT","WT_4_KID_GT","WT_5_KID_GT","WT_6_KID_GT")

vcf$CTX_GT_0.0 <- rowSums(vcf[,ctx]=="0/0")
vcf$HIP_GT_0.0 <- rowSums(vcf[,hip]=="0/0")
vcf$KID_GT_0.0 <- rowSums(vcf[,kid]=="0/0")


### Count number of samples with at least 1 supp. reads in the 0/0 samples of each tissue
sup_reads <- 1

vcf$CTX_THRESH_ALT_0.0 <- unlist( bplapply(1:nrow(vcf),getNumbReads,df=vcf,cols=ctx) )
vcf$HIP_THRESH_ALT_0.0 <- unlist( bplapply(1:nrow(vcf),getNumbReads,df=vcf,cols=hip) )
vcf$KID_THRESH_ALT_0.0 <- unlist( bplapply(1:nrow(vcf),getNumbReads,df=vcf,cols=kid) )



### write in bed format
chrom <- vcf$CHROM
start <- vcf$POS
end <- vcf$END
vcf$CHROM <- NULL
vcf$POS <- NULL
vcf$END <- NULL
vcf$SUPP_VEC <- NULL

vcf.final <- cbind( chrom, start, end, vcf )


### Handle duplicates (n=16)
### Among the duplicated SVs, select the ones with the highest mean of alternative reads
# Mark dupl
nrow1 <- nrow(vcf.final)
vcf.final$ID <- paste0(vcf.final$chrom,"-",vcf.final$start,"-",vcf.final$end,"-",vcf.final$SVTYPE)
dupl.id <- vcf.final[ duplicated(vcf.final$ID),"ID"]
vcf.final$is.dupl <- vcf.final$ID %in% dupl.id
# Get dupl. SVs
dupl <- vcf.final[vcf.final$is.dupl,]
dupl$mean_sup <- apply(dupl[, grepl("_ALT",colnames(dupl)) & !grepl("_ALT_",colnames(dupl)) ],1,mean)
# Select the lucky ones
dupl_selected <- as.data.frame(dupl %>% 
                                 group_by(ID) %>% 
                                 slice_max(mean_sup, n = 1,with_ties=F) )
dupl_selected$mean_sup <- NULL
# Remove duplicated from original file
n1 <- length(unique(vcf.final$ID))
vcf.final <- subset(vcf.final,is.dupl==FALSE)
# Add the selected duplicated
vcf.final <- rbind(vcf.final,dupl_selected)
rownames(vcf.final) <- 1:nrow(vcf.final)
# Sort
vcf.final <- vcf.final[order(vcf.final$chrom,vcf.final$start),]
# Check
n2 <- length(unique(vcf.final$ID))
if (n1!=n2) { stop("Some SVs have been lost") }
nrow2 <- nrow(vcf.final)
cat("\tdiscarded",nrow1-nrow2,"duplicated variants\n")

# Remove SV on scaffolds
vcf.final <- vcf.final[grepl("chr",vcf.final$chr),]
# discard scaffolds also from breakpoint2 in BND
vcf.final2 <- vcf.final[ is.na(vcf.final$breakpoint2) | grepl("chr",vcf.final$breakpoint2),]
# remove also chrM
vcf.final2 <- vcf.final2[!grepl("chrM",vcf.final2$breakpoint2),]

n3 <- nrow(vcf.final2)
rm(vcf.final)
cat("\tdiscarded",nrow2-n3,"variants on scaffolds/chrM\n")

#########################################################################################################


cat("Annotating SVs - genomic features\n")

# Read and format header of the input files
gene_anno <- fread(gene_annoN,data.table = F,header = F,sep='\t')
names(gene_anno) <- c("chr","start","end","gene_id","tx_name","strand","gene_name")
gene_anno <- gene_anno[,c("tx_name","gene_name")]
genome <- read.table(genomeN)
names(genome) <- c("chrom","length")


# Generate txdb object for annotation of genomic features
txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(gtf,organism="Mus musculus",chrominfo = genome))
feature_Grl <- suppressMessages(extractFeatureFromTxDb(TxDb=txdb,upstreamTSS=3000,downstreamTSS=3000,downstreamTTS=300))


# Parse vcf file to get a sample-by-sample SV list
vcf.final.sample <- getSampleSV(vcf.final2,reads_thresh = 1)
vcf.final.sample$sample_comparison <- vcf.final.sample$Sample



### Annotate differently BND and other variants - BND have two breakpoints to be annotated

# Annotate real SVs - non BND
vcf.final.sample.ann_noBND <- annotateGenomicFeatures(df=subset(vcf.final.sample,SVTYPE!="BND"),
                                                      priorities = priorities,
                                                      feature_Grl,
                                                      quiet = T,
                                                      is.shuffle = F,
                                                      txdb_raw = txdb)
# Annotate real SVs - BND
vcf.final.sample.ann_BND <- annotateGenomicFeaturesBND2(df=subset(vcf.final.sample,SVTYPE=="BND"),
                                                        priorities = priorities,
                                                        feature_Grl,
                                                        quiet = T,
                                                        is.shuffle = F,
                                                        txdb_raw = txdb)

# bind BND and non BND
vcf.final.sample.ann <- rbind(vcf.final.sample.ann_noBND,vcf.final.sample.ann_BND)
vcf.final.sample.ann$id <- NULL
vcf.final.sample.ann$sample_comparison <- NULL


# Add annotation to original vcf/bed file
vcf.final.sample.ann$id <- paste0(vcf.final.sample.ann$chrom,"-",vcf.final.sample.ann$start,"-",vcf.final.sample.ann$end)
unique.vcf.final.sample.ann <- unique(vcf.final.sample.ann[,c("id","tx_name","tx_coord",
                                                              "distanceTSS","distanceTx","annotation")])
vcf.final2$id <- paste0(vcf.final2$chrom,"-",vcf.final2$start,"-",vcf.final2$end)
vcf.final2 <- merge(vcf.final2,unique.vcf.final.sample.ann,by="id",all=T,sort=F)
vcf.final2$id <- NULL


### add gene symbol - do it diversely for BND and not BND as BND have a double tx annotation tx1/tx2

### non BND
vcf.final2_noBND <- subset(vcf.final2,SVTYPE!="BND")
vcf.final3_noBND <- getGeneSymbol(vcf.final2_noBND,gene_anno)

### BND
vcf.final2_BND <- subset(vcf.final2,SVTYPE=="BND")
# breakpoint 1
vcf.final2_BND_1 <- vcf.final2_BND
vcf.final2_BND_1$tx_name <- str_split_fixed(vcf.final2_BND_1$tx_name,"/",2)[,1]
vcf.final3_BND_1 <- getGeneSymbol(vcf.final2_BND_1,gene_anno)
# breakpoint 2
vcf.final2_BND_2 <- vcf.final2_BND
vcf.final2_BND_2$tx_name <- str_split_fixed(vcf.final2_BND_2$tx_name,"/",2)[,2]
vcf.final3_BND_2 <- getGeneSymbol(vcf.final2_BND_2,gene_anno)
# merge breakpoint1 and breakpoint2
vcf.final3_BND <- mergeBND(vcf.final2_BND,vcf.final3_BND_1,vcf.final3_BND_2,vcf.final2)

### Bind BND and non BND
vcf.final3 <- rbind(vcf.final3_noBND,vcf.final3_BND)

# order row as original
rownames(vcf.final3) <- vcf.final3$ID
vcf.final3 <- vcf.final3[vcf.final2$ID,]
# order col as original
vcf.final3 <- vcf.final3[,c(colnames(vcf.final2),"gene_name")]

n4 <- nrow(vcf.final3)
if (n4!=n3) {stop("some SVs were lost")}
#########################################################################################################


cat("Annotating SVs - TE overlap\n")


# Parse vcf file to get a sample-by-sample SV list
vcf.final.sample <- getSampleSV(vcf.final3,reads_thresh = 1)
vcf.final.sample$ID <- paste(vcf.final.sample$chrom,vcf.final.sample$start,vcf.final.sample$end,vcf.final.sample$SVTYPE,sep="-")


#### Overlap with TEs coords

# Non-BND
nonbnd <- subset(vcf.final3,SVTYPE!="BND")
write.table(nonbnd,"nonbnd",row.names = F,col.names = T,sep='\t',quote=F)
nonbnd.sample <-subset(vcf.final.sample,SVTYPE!="BND")
out_ls <- countTEoverlap2(path_dataf = "nonbnd",
                          dataf_melt = nonbnd.sample,
                          rmsk = rmsk,
                          bedtools = bedtools,
                          i = 6492)
system("rm nonbnd")
nonbnd_res <- out_ls[[2]]

# BND 1st breakpoint
bnd1 <- subset(vcf.final3,SVTYPE=="BND")
write.table(bnd1,"bnd1",row.names = F,col.names = T,sep='\t',quote=F)
bnd1.sample <-subset(vcf.final.sample,SVTYPE=="BND")
out_ls <- countTEoverlap2(path_dataf = "bnd1",
                          dataf_melt = bnd1.sample,
                          rmsk = rmsk,
                          bedtools = bedtools,
                          i = 6492)
system("rm bnd1")
bnd1_res <- out_ls[[2]]

# BND 2nd breakpoint
bnd2 <- bnd1
bnd2$chrom <- str_split_fixed(bnd2$breakpoint2,":",2)[,1]
bnd2$start <- as.numeric(str_split_fixed(bnd2$breakpoint2,":",2)[,2])
bnd2$end <- bnd2$start + 1
write.table(bnd2,"bnd2",row.names = F,col.names = T,sep='\t',quote=F)
out_ls <- countTEoverlap2(path_dataf = "bnd2",
                          dataf_melt = bnd1.sample,
                          rmsk = rmsk,
                          bedtools = bedtools,
                          i = 6492)
system("rm bnd2")
bnd2_res <- out_ls[[2]]

# merge BND1 and BND2
if (!all(bnd1_res$ID==bnd2_res$ID)) {stop("Some BNDs are found in BND1 but not BND2, or vice-versa")}
bnd <- merge(bnd1_res,bnd2_res[,c("ID","TE")],all=T,by="ID")
bnd$TE <- paste0(bnd$TE.x,"/",bnd$TE.y)
bnd$TE.x <- NULL
bnd$TE.y <- NULL



# Add annotation to original file
vcf.final.sample.ann <- rbind(nonbnd_res,bnd)
vcf.final.sample.ann$id <- paste0(vcf.final.sample.ann$chrom,"-",
                                  vcf.final.sample.ann$start,"-",vcf.final.sample.ann$end)
unique.vcf.final.sample.ann <- unique(vcf.final.sample.ann[,c("id","TE")])
vcf.final3$id <- paste0(vcf.final3$chrom,"-",
                        vcf.final3$start,"-",vcf.final3$end)
vcf.final4 <- merge(vcf.final3,unique.vcf.final.sample.ann,by="id",all=T,sort=F)
vcf.final4$id <- NULL

n5 <- nrow(vcf.final4)

if (n5!=n4) {stop("some SVs were lost when annotating over repeats")}
#########################################################################################################
# Overlap with simple repeats from tandem repeat finder (simple_rep1) and rmsk (simple_rep2)

cat("Annotating SVs - Simple repeats overlap\n")


### Non bnd
write.table(subset(vcf.final4,SVTYPE!="BND")[,c("chrom","start","end","ID")],"bedA.bed.tmp",row.names = F,col.names = T,sep='\t',quote=F)
res1 <- wrap_countSRoverlap(file = "bedA.bed.tmp",sr1 = simple_rep1,sr2 = simple_rep2,bedtools,i=6492)
intersect_SR1_nobnd <- res1[[1]] # overlap with TRF SR
names(intersect_SR1_nobnd)[ncol(intersect_SR1_nobnd)] <- "SR_ID1"
intersect_SR2_nobnd <- res1[[2]] # overlap with rmsk SR
names(intersect_SR2_nobnd)[ncol(intersect_SR2_nobnd)] <- "SR_ID2"


### BND
# BND breakpoint 1
write.table(subset(vcf.final4,SVTYPE=="BND")[,c("chrom","start","end","ID")],"bedA.bed.tmp",row.names = F,col.names = T,sep='\t',quote=F)
res2 <- wrap_countSRoverlap("bedA.bed.tmp",simple_rep1,simple_rep2,bedtools,i=6492)
intersect_SR1_bnd1 <- res2[[1]] # overlap with TRF SR
intersect_SR2_bnd1 <- res2[[2]] # overlap with rmsk SR

# BND breakpoint 2
nn <- subset(vcf.final4,SVTYPE=="BND")
nn$chrom <-str_split_fixed(nn$breakpoint2,":",2)[,1]
nn$start <- as.numeric(str_split_fixed(nn$breakpoint2,":",2)[,2])
nn$end <- nn$start+1
write.table(nn[,c("chrom","start","end","ID")],"bedA.bed.tmp",row.names = F,col.names = T,sep='\t',quote=F)
res3 <- wrap_countSRoverlap("bedA.bed.tmp",simple_rep1,simple_rep2,bedtools,i=6492)
intersect_SR1_bnd2 <- res3[[1]] # overlap with TRF SR
intersect_SR2_bnd2 <- res3[[2]] # overlap with rmsk SR

### Merge BND breakpoint1 and breakpoint2
# TRF SR
intersect_SR1_bnd <- mergeBND_SR(sr1_res = intersect_SR1_bnd1,sr2_res = intersect_SR1_bnd2,bp = 1)
# rmsk SR
intersect_SR2_bnd <- mergeBND_SR(sr1_res = intersect_SR2_bnd1,sr2_res = intersect_SR2_bnd2,2)

# bind BND and non-BND
intersect_SR1 <- rbind(intersect_SR1_nobnd,intersect_SR1_bnd) # TRF SR
intersect_SR2 <- rbind(intersect_SR2_nobnd,intersect_SR2_bnd) # rmsk SR

################################################################################################################################

# merge TRF SR and rmsk SR
intersect_SR <- merge(intersect_SR1,intersect_SR2,by="ID",all=T)

# Define whether selecting first or second SR annotation - priority to TRF
intersect_SR$SR <- ifelse(is.na(intersect_SR$SR_1),intersect_SR$SR_2, # if first is NA, get second
                          ifelse(is.na(intersect_SR$SR_2),intersect_SR$SR_1, # if second is NA, get first
                                 ifelse(intersect_SR$SR_1!=intersect_SR$SR_2,intersect_SR$SR_1, # if first is different from second, get first
                                        ifelse(intersect_SR$SR_1==intersect_SR$SR_2,intersect_SR$SR_1,"ERROR"))) )
# The same for intersect value
intersect_SR$intersect <- ifelse(is.na(intersect_SR$SR_1),intersect_SR$intersect_2, # if first is NA, get second
                          ifelse(is.na(intersect_SR$SR_2),intersect_SR$intersect_1, # if second is NA, get first
                                 ifelse(intersect_SR$SR_1!=intersect_SR$SR_2,intersect_SR$intersect_1, # if first is different from second, get first
                                        ifelse(intersect_SR$SR_1==intersect_SR$SR_2,intersect_SR$intersect_1,"ERROR"))) )

# The same for SR coord
intersect_SR$SR_ID <- ifelse(is.na(intersect_SR$SR_1),intersect_SR$SR_ID2, # if first is NA, get second
                                 ifelse(is.na(intersect_SR$SR_2),intersect_SR$SR_ID1, # if second is NA, get first
                                        ifelse(intersect_SR$SR_1!=intersect_SR$SR_ID2,intersect_SR$SR_ID1, # if first is different from second, get first
                                               ifelse(intersect_SR$SR_1==intersect_SR$SR_ID2,intersect_SR$SR_ID1,"ERROR"))) )
intersect_SR_f <- intersect_SR[,c("ID","SR","intersect","SR_ID")]
names(intersect_SR_f)[3] <- "intersect_SR"
intersect_SR_f$intersect_SR <- as.numeric(intersect_SR_f$intersect_SR)


# define which SR
n6 <- nrow(vcf.final4)
vcf.final5 <- merge(vcf.final4,intersect_SR_f,by="ID",all.x=T,sort=F)
n7 <- nrow(vcf.final5)
if (n6!=n7) {stop("some SVs were lost when annotating over simple repeats")}

# Calc intersection in %
vcf.final5$intersect_SR_perc <- vcf.final5$intersect_SR / vcf.final5$SVLEN * 100
vcf.final5$intersect_SR_perc <- round(vcf.final5$intersect_SR_perc,2)

### Adjust overlap with BND and INS

# if SVs are BND, set intersect_SR_perc to 100%
vcf.final5$intersect_SR_perc <- ifelse( (!is.na(vcf.final5$SR)) & vcf.final5$SVTYPE=="BND",100,vcf.final5$intersect_SR_perc)

# if SVs are INS, set intersect_SR_perc to 100% and intersect_SR to 1
vcf.final5$intersect_SR <- ifelse( (!is.na(vcf.final5$SR)) & vcf.final5$SVTYPE=="INS",1,vcf.final5$intersect_SR)
vcf.final5$intersect_SR_perc <- ifelse( (!is.na(vcf.final5$SR)) & vcf.final5$SVTYPE=="INS",100,vcf.final5$intersect_SR_perc)


# classify overlaps.simplerep 
if (table(intersect_SR_f$ID %in% vcf.final5$ID)["TRUE"]!=length(unique(intersect_SR_f$ID))) {
  stop("Some of the SVs in the intersect_SR_f are missing from the tot file")
}
vcf.final5$overlaps.simplerep <- vcf.final5$ID %in% intersect_SR_f$ID
#########################################################################################################


### Final annotation
vcf.final5$final.anno <- ifelse(vcf.final5$overlaps.simplerep==TRUE,"simple_rep","normal")

# Sort
vcf.final5 <- vcf.final5[order(vcf.final5$chrom,vcf.final5$start),]


#########################################################################################################

cat("Writing results to parsed_vcf.bed\n")

# Add ID as last col
ids <- vcf.final5$ID
vcf.final5$ID <- NULL
vcf.final5$ID <- ids
fwrite(vcf.final5,"parsed_vcf.bed",quote=F,sep='\t',row.names = F,na = NA)


cat("Done.\n")
cat("#############################################################\n")

