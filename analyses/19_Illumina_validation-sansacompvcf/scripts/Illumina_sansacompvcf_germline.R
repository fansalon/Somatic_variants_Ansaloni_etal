###########################################################################################################
### Validate sniffels2 SV with Illumina data - using sansa - germline
###########################################################################################################
rm(list = ls())
source("../functions/functions_nanopore.R")
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
###########################################################################################################

# vcf files
ont_vcf      <- "../../zenodo/input_data/WT_3_CTX.vcf"
ill_vcf_germ <- "../16_runManta_Illumina/output/germline/results/variants/diploidSV.vcf"
inp_ls <- list(ont_vcf, ill_vcf_germ)

cat(paste0("Reading vcf files:\n\t- ",ont_vcf,"\n\t- ",ill_vcf_germ,"\n"))

# bed anno ont
ont_vcf_anno  <- fread("../18_WT_3_CTX_parse_vcf/WT_3_CTX.vcf.bed",data.table=F,sep='\t',header=T)
ont_vcf_anno <- subset(ont_vcf_anno,MOSAIC!='MOSAIC')
colnames(ont_vcf_anno) <- paste0(colnames(ont_vcf_anno),"_ont")
colnames(ont_vcf_anno)[colnames(ont_vcf_anno)=="ID_ont"] <- "ID"
cat(paste0("Selecting ",nrow(ont_vcf_anno)," germline ONT\n"))

# bed anno illumina
ill_vcf_anno_germ  <- fread("../17_Illumina_parse_vcf/diploidSV.vcf.bed",data.table=F,sep='\t',header=T)
colnames(ill_vcf_anno_germ) <- paste0(colnames(ill_vcf_anno_germ),"_illumina")
colnames(ill_vcf_anno_germ)[colnames(ill_vcf_anno_germ)=="ID_illumina"] <- "ID"

###########################################################################################################


cat("Preparing input vcf files...\n")

for (file in inp_ls) {
  
  # cp locally
  system(paste0("cp ",file, " ."))
  
  # name
  filen <- basename(file)
  
  # bgzip
  system(paste0("bgzip ",filen))
  
  # generate index
  system(paste0("  bcftools index -t ",filen,".gz"))
  
}
###########################################################################################################

###########################################################################################################


cat(paste0("Comparing ",ont_vcf," and ",ill_vcf_germ,"\n"))

# get file names
ont_vcfN      <- basename(ont_vcf)
ill_vcf_germN <- basename(ill_vcf_germ)

# run sansa
system(paste0("sansa compvcf -m 0 -e 0 -o ont_vs_diploid -a ",ill_vcf_germN,".gz ",ont_vcfN,".gz"))

# read-in results
germl_ont <- fread("ont_vs_diploid.sv.classification",data.table=F,sep='\t',header=T)
germl_ont$BestMatchId <- gsub("base","",germl_ont$BestMatchId)

# annotate (somatic, simple repeats, TEs)
germl_ont_anno <- merge(germl_ont,ont_vcf_anno[,c("ID","SV_ID_ont","MOSAIC_ont","REF_reads_ont","ALT_reads_ont","MAF_ont",
                                                "overlaps.simplerep_ont","overlaps.TEs_ont")],by="ID",sort=F)
germl_ont_anno <- merge(germl_ont_anno,ill_vcf_anno_germ[,c("ID","SV_ID_illumina")],by.x="BestMatchId",by.y="ID",all.x=T,sort=F)

# sanity checks
if (!is.na(table(is.na(germl_ont_anno$SV_ID_ont))["TRUE"])) {stop("NA in SV_ID_ont")}
if (!table(is.na(germl_ont_anno$SV_ID_illumina))["TRUE"]==table(germl_ont_anno$Classification)["FP"]) {stop("NA in SV_ID_illumina")}

# final file
germl_ff <- germl_ont_anno[,c("ID","Classification","SVType","SV_ID_ont","SV_ID_illumina","MOSAIC_ont",
                              "REF_reads_ont","ALT_reads_ont","MAF_ont",
                              "overlaps.simplerep_ont","overlaps.TEs_ont")]

###########################################################################################################

###########################################################################################################

# clean
system("mkdir -p tmp/diploid/")
system("mv ont_vs_diploid* tmp/diploid")
system("mv diploid* tmp/diploid")
system("rm  WT*")
###########################################################################################################

###########################################################################################################


# write to out
fwrite(germl_ff,"ONT_germline_validation.txt",row.names = F,sep='\t',quote=F,na = NA)

# write summary to out
print(table(germl_ff[,c("Classification","SVType")]))
write.table(table(germl_ff[,c("Classification","SVType")]),"summary_germline.txt",sep='\t',quote=F)

