###########################################################################################################
### 10.
### Shuffling analysis - common somatic DELs - overlap with SR (N)n
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
###########################################################################################################


cat("Prepping files...\n")


#### Prepare file containing coords of the tissue specific private BNDs


# Select common somatic DELs in SR
common <- subset(bed,somatic.anno=="Somatic_common_private" | somatic.anno=="Somatic_common_recurrent" )
common <- subset(common,SVTYPE=="DEL" & overlaps.simplerep==T)

# Modify names in a nicer way ...
common$somatic.anno <- ifelse(common$somatic.anno=="Somatic_common_private","Somatic private",
                                               ifelse(common$somatic.anno=="Somatic_common_recurrent","Somatic recurrent","error"))

# get only (N)n
common$SR_motif_len <- nchar(common$SR) - 3 # do not count (, ) and n
common_2n <- subset(common,SR_motif_len==1)

###########################################################################################################


###########################################################################################################

### Get motif type (NN)n for the DELs of interest
obs <- getCountsBND_SR(dataf=common_2n,group="somatic.anno",item_to_count="SR")



############################################################################################################

### Get motifs type (NN)n for shuffle DELs


### Prep file to shuffle
# get new ID
common_2n$ID2 <- paste0(common_2n$somatic.anno,"@",common_2n$ID)
write.table(common_2n[,c("chrom","start","end","ID2")],"uno.tmp",col.names = F,row.names = F,sep='\t',quote=F)


### Prep file where to shuffle - SR (N)n - unique file
# get single SR file to shuffle the common DELs in - this must contain only SR (N)n
system(paste0("cat ",simple_rep1," ",simple_rep2," | sort | uniq | bedtools sort -i stdin > sr.bed"))
srbed <- fread("sr.bed",data.table=F)
srbed$len <- nchar(srbed$V4) -3
srbed <- subset(srbed,len==1)
srbed$len <- NULL
fwrite(srbed,"sr.bed",na = NA,row.names = F,col.names = F,sep='\t',quote=F)


### Prep individual SR file, get only (NN)n SR
simple_rep1_N <- getN_SR(simple_rep1,1)
simple_rep2_N <- getN_SR(simple_rep2,2)


# run function in parallel
cat("Started shuffling analysis...\n")
t0=Sys.time()
mergedls <- ( bplapply(seq_len(bootstrap),
                       BPPARAM=MulticoreParam(worker=multicoreWorkers(),RNGseed=seed_f),
                       getShuffleCommonDEL_SR_counts2_SR2,
                       common_2n,
                       genome,
                       "sr.bed",
                       simple_rep1_N,
                       simple_rep2_N,
                       bedtools,
                       bootstrap))
system("rm sr.bed uno.tmp sr1.bed.tmp sr2.bed.tmp")
t1=Sys.time()
t1-t0
# bootstrap=2       | 7''
# bootstrap=10      | 18''
# bootstrap=100     | 18''
# bootstrap=1000    | 1h

cat("Getting final results...\n")
res <- bindShuf_motif_common(inp_ls = mergedls,obs_df = obs,name="SR")
res$bootstrap <- bootstrap

write.table(res,paste0("DEL_somatic_common_SRmotif_shuffling_",bootstrap,".txt"),quote=F,sep='\t',row.names=F)


cat("Done\n")
cat("\n##################################################################################\n")

