###########################################################################################################
### 03.
### This script gets the SV from the 01 folder and the somatic annotation from 02 and annotate the SVs
### as common recurrent, common privare, tissue-specific recurrent or tissue-specific private
###########################################################################################################
rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(stringr))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################

cat("\n##################################################################################\n")
cat("Loading:\n\tvcf.bed file (parsed_vcf.bed)\n\tsomatic annotation object (somatic.Rds)\n")

# Load vcf
vcfbed <- fread("../01_parse_vcf/parsed_vcf.bed.gz",header=T,sep='\t',data.table = F)


# Load somatic SV
somatic_anno <- readRDS(("../02_get_somatic/somatic.Rds"))
somatic_common <- somatic_anno$common_somatic
somatic_tissue <- somatic_anno$tissue_specific_somatic


cat("Getting ID of:\n\tcommon private SVs\n\tcommon recurrent SVs\n\ttissue-specific private SVs\n\ttissue-specific recurrent SVs\n")


### Define if somatic SVs are private or recurrent

### Common
# count how many tissues present each SV in each mouse
somatic_common$n_samples <- apply(somatic_common,1,sum)
# define private or recurrent
somatic_common_private <- rownames(subset(somatic_common,n_samples==1))
#6861
##somatic_common_recurrent <- rownames(subset(somatic_common,n_samples>1))
##26790
somatic_common_recurrent <- rownames(subset(somatic_common,n_samples==5))
#8164
somatic_common_unclassified <- rownames(subset(somatic_common,n_samples>1 & n_samples<5))
#18612


### Tissue-specific
# count how many tissues present each SV in each mouse
somatic_tissue$n_samples <- apply(somatic_tissue,1,sum)
# define private or recurrent
somatic_tissue_private <- rownames(subset(somatic_tissue,n_samples==1))
#77615
##somatic_tissue_recurrent <- rownames(subset(somatic_tissue,n_samples>1))
###4994
somatic_tissue_recurrent <- rownames(subset(somatic_tissue,n_samples==5))
#325
somatic_tissue_unclassified <- rownames(subset(somatic_tissue,n_samples>1 & n_samples<5))
#4701



# Annotate
cat("Annotating original SVs\n")
vcfbed$somatic.anno <- ifelse(vcfbed$ID %in% somatic_common_private,"Somatic_common_private",
                              ifelse(vcfbed$ID %in% somatic_common_recurrent,"Somatic_common_recurrent",
                                     ifelse(vcfbed$ID %in% somatic_common_unclassified,"Somatic_common_unclassified",
                                            ifelse(vcfbed$ID %in% somatic_tissue_private,"Somatic_tissue_private",
                                                   ifelse(vcfbed$ID %in% somatic_tissue_recurrent,"Somatic_tissue_recurrent",
                                                          ifelse(vcfbed$ID %in% somatic_tissue_unclassified,"Somatic_tissue_unclassified",NA))))))
# checks
if (nrow(subset(vcfbed,somatic.anno=="Somatic_common_private"))!=length(somatic_common_private)){stop()}
if (nrow(subset(vcfbed,somatic.anno=="Somatic_common_recurrent"))!=length(somatic_common_recurrent)){stop()}
if (nrow(subset(vcfbed,somatic.anno=="Somatic_common_unclassified"))!=length(somatic_common_unclassified)){stop()}
if (nrow(subset(vcfbed,somatic.anno=="Somatic_tissue_private"))!=length(somatic_tissue_private)){stop()}
if (nrow(subset(vcfbed,somatic.anno=="Somatic_tissue_recurrent"))!=length(somatic_tissue_recurrent)){stop()}
if (nrow(subset(vcfbed,somatic.anno=="Somatic_tissue_unclassified"))!=length(somatic_tissue_unclassified)){stop()}
#############################################################################################################################


# Define samples where SV is found and concat this into one string
cat("Defining sample(s) where each SV is found\n")

# Common


for (i in 1:nrow(vcfbed)) {
  
  # counter
  if (i%%10000==0) { cat(paste0("\t",i," rows parsed (",round(i/nrow(vcfbed)*100,0),"%)\n")) }
  
  # get SV row by row
  d1 <- vcfbed[i,]
  
  # if somatic!
  if (!is.na(d1$somatic.anno)) {
    
    
    if (d1$somatic.anno=="Somatic_common_recurrent") {
      
      # count in how many tissue the SV is found in each mouse
      matches2_tx <- getsample(d1)
      
      # these are common recurrent SVs, thus are found in 3/3 tissues...
      samples <- names(matches2_tx[matches2_tx==3])
      # ..in 5 mice!
      if (length(samples)!=5) {stop("These are common recurrent SVs")}
      
      
    } else if (d1$somatic.anno=="Somatic_common_private") {
      
      # count in how many tissue the SV is found in each mouse
      matches2_tx <- getsample(d1)
      
      # these are common private SVs, thus are found in 3/3 tissues..
      samples <- names(matches2_tx[matches2_tx==3])
      # ..in 1 mouse!
      if (length(samples)!=1) {stop("These are common private SVs")}
      
      
    } else if (d1$somatic.anno=="Somatic_common_unclassified") {
      
      # count in how many tissue the SV is found in each mouse
      matches2_tx <- getsample(d1)
      
      # these are common unclassified SVs, thus are found in 3/3 tissues..
      samples <- names(matches2_tx[matches2_tx==3])
      # ..in 2-to-4 mice
      if (length(samples)==1 | length(samples)==5) {stop("These are common private SVs")}
      
      
    } else if (d1$somatic.anno=="Somatic_tissue_private") {
      
      # these are tissue-specific private SVs, thus must be found in only 1 sample
      samples <- getsampleTx(d1)
      mice <- length(unique(paste0("WT_",str_split_fixed(samples,"_",3)[,2])))
      if (mice!=1) {
        cat("These are tissue-specific private SVs\n1 mouse expected.\nGot:\n")
        print(samples)
        stop("exiting")
      }
      # modify name
      samples <- paste0(str_split_fixed(samples,"_",3)[,3],"_",str_split_fixed(samples,"_",3)[,2])
      
      
    } else if (d1$somatic.anno=="Somatic_tissue_recurrent") {
      
      # these are tissue-specific recurrent SVs, thus must be found in only 1 tissue, but >=5 mice
      samples <- getsampleTx(d1)
      if (length(samples)<5) {stop("These are tissue-specific recurrent SVs")}
      # modify name
      samples <- paste0(str_split_fixed(samples,"_",3)[,3],"_",str_split_fixed(samples,"_",3)[,2])
      
      
    } else if (d1$somatic.anno=="Somatic_tissue_unclassified") {
      
      # these are tissue-specific unclassified SVs, thus must be found in only 1 sample, in at least two mouse
      samples <- getsampleTx(d1)
      if (length(samples)<=1) {stop("These are tissue-specific unclassified SVs")}
      # modify name
      samples <- paste0(str_split_fixed(samples,"_",3)[,3],"_",str_split_fixed(samples,"_",3)[,2])
      
      
    } else {
      cat("SV not identified in any known class\n")
      print(d1)
      stop("exiting")
    }
    # convert character to 
    vcfbed[i,"Samples"] <- paste(rev(samples),sep="",collapse=",")
    
    
  }
}
#############################################################################################################################



cat("Writing output file (parsed_vcf_anno_somatic.bed)\n")
fwrite(vcfbed,"parsed_vcf_anno_somatic.bed",row.names = F,sep='\t',quote=F,na = NA)

cat("Done.\n")
cat("##################################################################################\n")

