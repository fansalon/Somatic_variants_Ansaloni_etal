###########################################################################################################
### 07.
### Shuffling analysis - tissue specific recurrent somatic BND - overlap with SR (NN)n
###########################################################################################################
rm(list = ls())
seed_f='110322'
bedtools="/opt/homebrew/bin/bedtools"
genome="../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr"
simple_rep1="../../zenodo/input_data/simpleRepeat.txt.bed"
simple_rep2="../../zenodo/input_data/simpleRepeat.rmsk.bed"
bootstrap=1000
###########################################################################################################
source("../functions/functions_nanopore.R")
library(BiocParallel)
library(data.table)
library(dplyr)
library(stringr)
library(rlist)
###########################################################################################################

cat("\n##################################################################################\n")
cat("Loading:\n\t../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz\n")

# Somatic anno+features
bed <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)

# Load somatic obj
cat("Loading:\n\t../02_get_somatic/somatic.Rds\n")
obj <- readRDS("../02_get_somatic/somatic.Rds")
comparison <- obj$tissue_specificity
comparison$Type <- gsub("_vs_","-vs-",comparison$Type)

###########################################################################################################


cat("Prepping files...\n")


#### Prepare file containing coords of the tissue specific private BNDs


# Select tissue-specific somatic BNDs in SR
tissue_specific <- subset(bed,somatic.anno=="Somatic_tissue_private" & SVTYPE=="BND" & overlaps.simplerep==T)


# Add pairwise comparison type
tissue_specific2 <- merge(tissue_specific,comparison,by.x="ID",by.y="SV_ID",all.x=T,sort=F)

# Expand dataframe based on the "Samples" column (which contains a comma-separated list of the samples where the sv is found)
# if a single SV is associated to n sample, report it n times, each time associated to a different sample
tissue_specific2_expand <- as.data.frame(tidyr::separate_rows(tissue_specific2, Samples, sep=","))

# Discard samples not matching comparison type
# there might be cases where the samples are KID_6,CTX_2 and the comparisons KID-vs-HIP and CTX-vs-HIP
# in such cases, when the df is expanded, we need to select the matching samples and comparisons
tissue_specific2_expand$tissue_sample <- str_split_fixed(tissue_specific2_expand$Samples,"_",2)[,1]
tissue_specific2_expand$tissue_type <- str_split_fixed(tissue_specific2_expand$Type,"-",3)[,1]
tissue_specific2_expand <- subset(tissue_specific2_expand,tissue_sample==tissue_type)

# Define non-redundant comparison type
tissue_specific2_expand$Type2 <- 
  ifelse(tissue_specific2_expand$Type=="HIP-vs-KID" |tissue_specific2_expand$Type=="KID-vs-HIP","HIP vs\nKID",
         ifelse(tissue_specific2_expand$Type=="CTX-vs-KID" |tissue_specific2_expand$Type=="KID-vs-CTX","CTX vs\nKID",
                ifelse(tissue_specific2_expand$Type=="HIP-vs-CTX" |tissue_specific2_expand$Type=="CTX-vs-HIP","CTX vs\nHIP","error")))

# Modify names in a nicer way ...
tissue_specific2_expand$somatic.anno <- ifelse(tissue_specific2_expand$somatic.anno=="Somatic_tissue_private","Somatic private",
                                               ifelse(tissue_specific2_expand$somatic.anno=="Somatic_tissue_recurrent","Somatic recurrent","error"))

# These are BND, there are two breakpoints that can be both, or only one, associated to SR, separate them!
# BP 1
rep_1 <- tissue_specific2_expand
rep_1$SR <- str_split_fixed(rep_1$SR,"/",2)[,1]
rep_1 <- subset(rep_1,SR!="No SR")
rep_1 <- rep_1[,c("chrom","start","end","Type","ID","SR")]
rep_1$bp <- "bp_1"

# BP 2
rep_2 <- tissue_specific2_expand
rep_2$SR <- str_split_fixed(rep_2$SR,"/",2)[,2]
rep_2 <- subset(rep_2,SR!="No SR")
rep_2 <- rep_2[,c("chrom","start","end","Type","ID","breakpoint2","SR")]
rep_2$bp <- "bp_2"
rep_2$chrom <- str_split_fixed(rep_2$breakpoint2,":",2)[,1]
rep_2$start <- as.numeric(str_split_fixed(rep_2$breakpoint2,":",2)[,2])
rep_2$end <- rep_2$start + 1
rep_2$breakpoint2 <- NULL


rep <- rbind(rep_1,rep_2)
# get only (NN)n
rep$SR_motif_len <- nchar(rep$SR) - 3 # do not count (, ) and n
rep <- subset(rep,SR_motif_len==2)

###########################################################################################################


###########################################################################################################

### Get motif type (NN)n for the BND of intetest
obs <- getCountsBND_SR(dataf=rep,group="Type",item_to_count="SR")



############################################################################################################

### Get motifs type (NN)n  for shuffle BNDs


### Prep file to shuffle
# get new ID
rep$ID2 <- paste0(rep$Type,"@",rep$ID,"@",rep$bp)
write.table(rep[,c("chrom","start","end","ID2")],"uno.tmp",col.names = F,row.names = F,sep='\t',quote=F)


### Prep file where to shuffle - SR (NN)n - unique file
# get single SR file to shuffle the BND in - this must contain only SR (NN)n
system(paste0("cat ",simple_rep1," ",simple_rep2," | sort | uniq | bedtools sort -i stdin > sr.bed"))
srbed <- fread("sr.bed",data.table=F)
srbed$len <- nchar(srbed$V4) -3
srbed <- subset(srbed,len==2)
srbed$len <- NULL
fwrite(srbed,"sr.bed",na = NA,row.names = F,col.names = F,sep='\t',quote=F)


### Prep individual SR file, get only (NN)n SR
simple_rep1_NN <- getNN_SR(simple_rep1,1)
simple_rep2_NN <- getNN_SR(simple_rep2,2)


# run function in parallel
cat("Started shuffling analysis...\n")
t0=Sys.time()
mergedls <- ( bplapply(seq_len(bootstrap),
                       BPPARAM=MulticoreParam(worker=multicoreWorkers(),RNGseed=seed_f),
                       getShuffleBND_SR_counts2_SR2,
                       rep,
                       genome,
                       "sr.bed",
                       simple_rep1_NN,
                       simple_rep2_NN,
                       bedtools,
                       bootstrap))
system("rm sr.bed uno.tmp sr1.bed.tmp sr2.bed.tmp")
t1=Sys.time()
t1-t0
# bootstrap=2       | 16''
# bootstrap=10      | 37''
# bootstrap=1000    | 1h

cat("Getting final results...\n")
res <- bindShuf_motif(inp_ls = mergedls,obs_df = obs,name="SR")
res$bootstrap <- bootstrap

write.table(res,paste0("BND_somatic_private_SRmotiff_shuffling_",bootstrap,".txt"),quote=F,sep='\t',row.names=F)


cat("Done\n")
cat("\n##################################################################################\n")
