###########################################################################################################
### 02.
### This script parses the vcf.bed file extracting somatic SVs - this prepares figure panels for Fig1 and Fig2
###########################################################################################################
rm(list = ls())
cat("\n##################################################################################\n")
cat("Loading libraries\n")
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(scales))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(ggridges))
#######################################################################################################################################
source("../functions/functions_nanopore.R")
#######################################################################################################################################
#######################################################################################################################################
cat("Loading vcf.bed file (parsed_vcf.bed)\n")
vcf.final <- fread("../01_parse_vcf/parsed_vcf.bed.gz",header=T,sep='\t',data.table = F)
rownames(vcf.final) <- vcf.final$ID
out_res <- list()
#######################################################################################################################################


#### Define SVs as 1/1, 0/1 and 0/0
## 1/1 is prioritized over 0/1 and 0/0 - this means that SV found as 1/1 in at least 1 sample are classified as 1/1, and so on

cat("Classifying the",nrow(vcf.final),"SVs as 0/1, 1/1 or 0/0\n")

# Get columns reporting genotype
genotype <- vcf.final[,grepl("_GT",colnames(vcf.final)) & !grepl("_0.0",colnames(vcf.final))]
nco <- ncol(genotype)
# count how many times the SV is genotyped 0/0, 0/1 or 1/1
genotype$GT_00 <- rowSums(genotype=="0/0")
genotype$GT_01 <- rowSums(genotype=="0/1")
genotype$GT_11 <- rowSums(genotype=="1/1")
genotype$GT_count <- ifelse(genotype$GT_11>=1,"1/1",
                            ifelse(genotype$GT_01>=1,"0/1",
                                   ifelse(genotype$GT_00==nco,"0/0","./.")))
# discard ./.
genotype <- subset(genotype,GT_count!="./.")
genotype_count <- as.data.frame(table(genotype$GT_count ))
genotype_count <- genotype_count[order(genotype_count$Var1,decreasing = T),]
genotype_count$Var1 <- factor(genotype_count$Var1,levels = unique(genotype_count$Var1))

### --> used for Fig 1C
out_res[["GT_count"]] <- genotype_count

#######################################################################################################################################


#######################################################################################################################################


#### Get somatic SVs only and define in how many tissues the SV is found (0/3, 1/3, 2/3, 3/3)

cat("Counting in how many tissues each SV is found in each individual mouse\n")

# get 0/0 --> SVs that are always 0/0, in all the 15 samples
somatic <- subset(genotype,GT_count=="0/0")

# get ALT reads of the 0/0
somatic_ALT <- vcf.final[rownames(somatic),grepl("_ALT",colnames(vcf.final)) & !grepl("THRESH",colnames(vcf.final))]

# count in how many tissues the SV is found (>=1 read) in each mouse
for (i in 2:6) {
  
  # Get ALT reads per each SV in each mouse (3 tissues)
  mouse <- somatic_ALT[,grepl(paste0("WT_",i,"_"),colnames(somatic_ALT))]
  # Count in how many tissues the SV is found
  mouse[,paste0("n_mouse_WT",i)] <- rowSums(mouse>=1)
  
  # Count how many SVs are found in 0,1,2,3 mice
  count <- as.data.frame(table(mouse[,paste0("n_mouse_WT",i)]))
  count$Mouse <- paste0("Mouse ",i-1)
  
  # Get dataf with SV and number of mice it is found in
  merged <- as.data.frame(cbind(mouse[,paste0("n_mouse_WT",i)],rownames(mouse)))
  names(merged) <- c(paste0("n_mouse_WT",i),"SV")
  merged[,1] <- as.numeric(merged[,1])
  
  # bind / merge
  if (i==2) {
    final_count <- count
    ff <- merged
  } else {
    final_count <- rbind(final_count,count)
    ff <- merge(ff,merged,by="SV",all=T)
  }
}
# discard SV never found 
final_count <- subset(final_count,Var1!="0")

# count total somatic SV per mouse
totals <- aggregate(final_count$Freq,by=list(final_count$Mouse),sum)
names(totals) <- c("Mouse","Tot")
final_count.aggr <- merge(final_count,totals,by="Mouse",all=T,sort=F)
final_count.aggr$Var1 <- paste0(final_count.aggr$Var1,"/3")

### --> used for Fig 1D
out_res[["somaticSV_mouse_tissue_count"]] <- final_count.aggr

#######################################################################################################################################


#######################################################################################################################################


#### Get common (3/3 tissues in at least 1 mouse) and tissue specific SVs


cat("Defining common and tissue-specific SVs\n")

rownames(ff) <- ff$SV
ff$SV <- NULL


#### Common
# count in how many mice the SVs is found in 3 tissues
ff$n_3.3 <- rowSums(ff==3)

# Get common (SVs found in 3/3 tissues in at least 1 mouse)
common_somatic <- subset(ff,n_3.3>=1)
common_somatic$n_3.3 <- NULL

# convert to binary (1 if 3/3, 0 if not)
for (j in 1:ncol(common_somatic)) {
  common_somatic[,j] <- ifelse(common_somatic[,j]==3,1,0)
}


#### Tissue-specific
# Get all the other SVs
tissue_somatic <- subset(ff,n_3.3<1)
tissue_somatic$n_3.3 <- NULL

# check the SVs are always found in at least one mouse
tissue_somatic2 <- tissue_somatic[apply(tissue_somatic,1,sum)>=1,]


################################################### *** IMPORTANT *** ###################################################
# tissue-specificity means that the SVs must not be found in all the 3 tissues, even among different mice               #
# count in how many tissues the SV is found                                                                             #
#########################################################################################################################

# Get ALT reads for the tissue-specific SVs
tissue_somatic2_ALT <- vcf.final[rownames(tissue_somatic2),grepl("_ALT",colnames(vcf.final)) & !grepl("THRE",colnames(vcf.final)) ]
# Set to TRUE if the SV is found in at least 1 ctx, hip and kid tissues
tissue_somatic2_ALT$CTX <- rowSums(tissue_somatic2_ALT[,grepl("CTX",colnames(tissue_somatic2_ALT))]>=1)>=1
tissue_somatic2_ALT$HIP <- rowSums(tissue_somatic2_ALT[,grepl("HIP",colnames(tissue_somatic2_ALT))]>=1)>=1
tissue_somatic2_ALT$KID <- rowSums(tissue_somatic2_ALT[,grepl("KID",colnames(tissue_somatic2_ALT))]>=1)>=1
# convert F/T to binary
tissue_somatic2_ALT$CTX <- ifelse(tissue_somatic2_ALT$CTX==T,1,0)
tissue_somatic2_ALT$HIP <- ifelse(tissue_somatic2_ALT$HIP==T,1,0)
tissue_somatic2_ALT$KID <- ifelse(tissue_somatic2_ALT$KID==T,1,0)
# count in how many tissues the SV is foumd
tissue_somatic2_ALT$n_tissues <- apply(tissue_somatic2_ALT[,c("CTX","KID","HIP")],1,sum)

# discard SVs found in 3/3 tissues - this are supposed to be tissue-specific SVs
# if a SV is found in ctx in mouse1, hip in mouse2 and kid in mouse3, it is not defenitely tissue-specific
# but it is not even common, as it is never found in all the three tissues of a single individual mouse
tissue_somatic_final <- subset(tissue_somatic2_ALT,n_tissues!=3)
tissue_somatic_final$n_tissues <- NULL
tissue_somatic_final$CTX <- NULL
tissue_somatic_final$HIP <- NULL
tissue_somatic_final$KID <- NULL
discarded <- nrow(tissue_somatic2_ALT)-nrow(tissue_somatic_final)

# convert to binary (1 if found, 0 not found)
for (j in 1:ncol(tissue_somatic_final)) {
  tissue_somatic_final[,j] <- ifelse(tissue_somatic_final[,j]>=1,1,0)
}

# Define tissue-specificity by making tissue pairwise comparisons
tissues <- c("CTX","HIP","KID")
res_list <- list()
for (i in 1:length(tissues)) {
  for (j in 1:length(tissues)) {
    if (i!=j) {
      # Get tissues
      tx1 <- tissues[i]
      tx2 <- tissues[j]
      comparison_name <- paste0(tx1,"_vs_",tx2)
      print(comparison_name)
      # Get ALT reads for all the mice, for only these tissues
      comp <- tissue_somatic_final[,grepl(tx1,colnames(tissue_somatic_final)) | grepl(tx2,colnames(tissue_somatic_final))]      # count how many time the SV is found in each tissue
      tx1_count <- rowSums(comp[,grepl(tx1,colnames(comp))>=1])
      tx2_count <- rowSums(comp[,grepl(tx2,colnames(comp))>=1])

      # get the tissue-specific
      # to be tissue-specific a given SV must be found in at least 1 sample in the first tissue
      # and never found in the second tissue
      binded <- cbind(as.data.frame(tx1_count),as.data.frame(tx2_count))
      tx_spec <- rownames(subset(binded,tx1_count>=1 & tx2_count==0))
      res_list[[comparison_name]] <- tx_spec
    }
  }
}


### Tissue-specific
# count how many tissues present each SV in each mouse (either 0, 1 or 2)
for (i in 2:6) {
  ind_mouse <- tissue_somatic_final[,grepl(i,colnames(tissue_somatic_final))]
  tissue_somatic_final[,paste0("Mouse_",i)] <- rowSums(ind_mouse>=1)
}
# convert to binary (1 if found in 1 or 2 tissues, 0 otherwise)
tissue_somatic_final2 <- tissue_somatic_final[,grepl("Mouse",colnames(tissue_somatic_final))]
for (j in 1:ncol(tissue_somatic_final2)) {
  tissue_somatic_final2[,j] <- ifelse(tissue_somatic_final2[,j]>0,1,0)
}


### --> used for Fig 2A and 2B
out_res[["common_somatic"]] <- common_somatic
out_res[["tissue_specific_somatic"]] <- tissue_somatic_final2
#######################################################################################################################################

#######################################################################################################################################

### Now create a new obj with the tissue-specific annotation
for (i in 1:length(res_list)) {
  d1 <- data.frame("SV_ID"=res_list[[i]],"Type"=names(res_list)[i])
  if (i==1) {
    tot <- d1
  } else {
    tot <- rbind(d1,tot)
  }
}
out_res[["tissue_specificity"]] <- tot



cat("Discarded",discarded,"tissue-specific SVs found in 3/3 tissues, in different mice\n")
cat("\tcommon SVs:",nrow(common_somatic),"\n")
cat("\ttissue-specific SVs:",nrow(tissue_somatic_final),"\n")

cat("Writing results to somatic.Rds\n")
saveRDS(out_res,"somatic.Rds")


cat("Done.\n")
cat("##################################################################################\n")

