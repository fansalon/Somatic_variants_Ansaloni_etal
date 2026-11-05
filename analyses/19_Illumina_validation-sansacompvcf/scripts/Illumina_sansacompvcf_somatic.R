###########################################################################################################
### Validate sniffels2 SV with Illumina data - using sansa - somatic
###########################################################################################################
###########################################################################################################
rm(list = ls())
source("../functions/functions_nanopore.R")
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
###########################################################################################################

# vcf files
ont_vcf      <- "../../zenodo/input_data/WT_3_CTX.vcf"
ill_vcf_somt <- "../16_runManta_Illumina/output/somatic/results/variants/tumorSV_mod.vcf"
inp_ls <- list(ont_vcf, ill_vcf_somt)

cat(paste0("Reading vcf files:\n\t- ",ont_vcf,"\n\t- ",ill_vcf_somt,"\n"))

# bed anno ont
ont_vcf_anno  <- fread("../18_WT_3_CTX_parse_vcf/WT_3_CTX.vcf.bed",data.table=F,sep='\t',header=T)
ont_vcf_anno <- subset(ont_vcf_anno,MOSAIC=='MOSAIC')
colnames(ont_vcf_anno) <- paste0(colnames(ont_vcf_anno),"_ont")
colnames(ont_vcf_anno)[colnames(ont_vcf_anno)=="ID_ont"] <- "ID"
cat(paste0("Selecting ",nrow(ont_vcf_anno)," somatic ONT\n"))

# bed anno illumina
ill_vcf_anno_smtc  <- fread("../17_Illumina_parse_vcf/tumorSV.vcf.bed",data.table=F,sep='\t',header=T)
colnames(ill_vcf_anno_smtc) <- paste0(colnames(ill_vcf_anno_smtc),"_illumina")
colnames(ill_vcf_anno_smtc)[colnames(ill_vcf_anno_smtc)=="ID_illumina"] <- "ID"

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


cat(paste0("Comparing ",ont_vcf," and ",ill_vcf_somt,"\n"))

# get file names
ont_vcfN      <- basename(ont_vcf)
ill_vcf_somtN <- basename(ill_vcf_somt)

# run sansa
system(paste0("sansa compvcf -m 0 -e 0 -o ont_vs_tumor -a ",ill_vcf_somtN,".gz ",ont_vcfN,".gz"))

# read-in results
smtc_ont <- fread("ont_vs_tumor.sv.classification",data.table=F,sep='\t',header=T)
smtc_ont$BestMatchId <- gsub("base","",smtc_ont$BestMatchId)

# annotate (somatic, simple repeats, TEs)
smtc_ont_anno <- merge(smtc_ont,ont_vcf_anno[,c("ID","SV_ID_ont","MOSAIC_ont","REF_reads_ont","ALT_reads_ont","MAF_ont",
                                                "overlaps.simplerep_ont","overlaps.TEs_ont")],by="ID",sort=F)
smtc_ont_anno <- merge(smtc_ont_anno,ill_vcf_anno_smtc[,c("ID","SV_ID_illumina")],by.x="BestMatchId",by.y="ID",all.x=T,sort=F)

# sanity checks
if (!is.na(table(is.na(smtc_ont_anno$SV_ID_ont))["TRUE"])) {stop("NA in SV_ID_ont")}
if (!table(is.na(smtc_ont_anno$SV_ID_illumina))["TRUE"]==table(smtc_ont_anno$Classification)["FP"]) {stop("NA in SV_ID_illumina")}

# final file
smtc_ff <- smtc_ont_anno[,c("ID","Classification","SVType","SV_ID_ont","SV_ID_illumina","MOSAIC_ont",
                            "REF_reads_ont","ALT_reads_ont","MAF_ont",
                            "overlaps.simplerep_ont","overlaps.TEs_ont")]

###########################################################################################################

###########################################################################################################


cat("Validating BND..\n")

# validate BND
validated_bnd <- validatingBND(ont=ont_vcf_anno, ill=ill_vcf_anno_smtc, D=1000)

# add validated to original file BND
for (i in 1:nrow(smtc_ff)) {
  if (smtc_ff[i,"SV_ID_ont"] %in% validated_bnd$SV_ID_ont) {
    val_ill_id <- validated_bnd[validated_bnd$SV_ID_ont==smtc_ff[i,"SV_ID_ont"],"SV_ID_illumina"]
    smtc_ff[i,"SV_ID_illumina"] <- val_ill_id
    smtc_ff[i,"Classification"] <- "TP"
  }
}

###########################################################################################################

###########################################################################################################

# clean
system("mkdir -p tmp/tumor/")
system("mv ont_vs_tumor* tmp/tumor")
system("mv tumor* tmp/tumor")
system("rm  WT*")
###########################################################################################################

###########################################################################################################


# write to out
fwrite(smtc_ff,"ONT_somatic_validation.txt",row.names = F,sep='\t',quote=F,na = NA)

# write summary to out
print(table(smtc_ff[,c("Classification","SVType")]))
write.table(table(smtc_ff[,c("Classification","SVType")]),"summary_somatic.txt",sep='\t',quote=F)

