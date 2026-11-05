###########################################################################################################
### 09.
### Shuffling analysis - common DELs - overlap with SR
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


#### Prepare file containing coords of the common DELs


# Select tissue-specific somatic BNDs in SR
common <- subset(bed,grepl("Somatic_comm",bed$somatic.anno) & SVTYPE=="DEL" & overlaps.simplerep==T)
common <- subset(common,somatic.anno!="Somatic_common_unclassified")

# Modify names in a nicer way ...
common$somatic.anno <- ifelse(common$somatic.anno=="Somatic_common_private","Somatic private",
                                               ifelse(common$somatic.anno=="Somatic_common_recurrent","Somatic recurrent","error"))

# add SR motif len
common$SR_motif_len <- nchar(common$SR) - 3 # -3 as we dont want to count "(", ")" and "n"
###########################################################################################################


###########################################################################################################

### Get motif lenght counts for the BND of intetest
obs <- getCountsBND_SR(dataf=common,group="somatic.anno",item_to_count="SR_motif_len")

############################################################################################################

### Get motif lenght counts for shuffle BNDs

### Prep file to shuffle
# get new ID
common$ID2 <- paste0(common$somatic.anno,"@",common$ID)
# write to file
write.table(common[,c("chrom","start","end","ID2")],"uno.tmp",col.names = F,row.names = F,sep='\t',quote=F)

### Prep where to shuffle
# get single SR file to shuffle the BND in
system(paste0("cat ",simple_rep1," ",simple_rep2," | sort | uniq | bedtools sort -i stdin > sr.bed"))


# run function in parallel
cat("Started shuffling analysis...\n")
t0=Sys.time()
mergedls <- ( bplapply(seq_len(bootstrap),
                       BPPARAM=MulticoreParam(worker=multicoreWorkers(),RNGseed=seed_f),
                       getShuffleBND_SR_counts2_common,
                       common,
                       genome,
                       "sr.bed",
                       simple_rep1,
                       simple_rep2,
                       bedtools,
                       bootstrap))
system("rm sr.bed uno.tmp")
t1=Sys.time()
t1-t0
# bootstrap=2       | 1'
# bootstrap=10      | 3'
# bootstrap=1000    | xh

cat("Getting final results...\n")
res <- bindShuf_common(mergedls,obs,"SR_motif_len")
res$id <- NULL
res$bootstrap <- bootstrap

write.table(res,paste0("DEL_somatic_common_SR_shuffling_",bootstrap,".txt"),quote=F,sep='\t',row.names=F)


cat("Done\n")
cat("\n##################################################################################\n")
