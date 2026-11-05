###########################################################################################################
### 05.
### GO analyses on variants of interest
###########################################################################################################
rm(list = ls())
###########################################################################################################
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(gprofiler2))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################
GP_VERSION <- "e111_eg58_p18"
gprofiler2::set_base_url(gp_archive_url(GP_VERSION))
vi <- tryCatch(gprofiler2::get_version_info("mmusculus"), error = function(e) NULL)
if (!is.null(vi)) cat("g:Profiler:", vi$gprofiler_version,
                      "|", vi$sources$`GO:BP`$version, "\n")
###########################################################################################################

cat("\n##################################################################################\n")
cat("Loading:\n\t../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz\n")

# Load vcf
vcfbed <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
# split annotations
vcfbed$annotation1 <- str_split_fixed(vcfbed$annotation,"/",2)[,1]
vcfbed$annotation2 <- str_split_fixed(vcfbed$annotation,"/",2)[,2]
vcfbed$gene_name1 <- str_split_fixed(vcfbed$gene_name,"/",2)[,1]
vcfbed$gene_name2 <- str_split_fixed(vcfbed$gene_name,"/",2)[,2]


# Load somatic obj
cat("Loading:\n\t../02_get_somatic/somatic.Rds\n")
obj <- readRDS("../02_get_somatic/somatic.Rds")
###########################################################################################################





cat("\n################################# Running GO on all the somatic variants")
# background is any gene hosting a somatic variant
# (within the function only genes hosting genic variants are selected)
somatic <- subset(vcfbed,!is.na(vcfbed$somatic.anno))
# run go
final_res_all <- wrapGO_tissue_specific(background=somatic,
                                        svtype="BND",
                                        svclass="Somatic_tissue_private",
                                        fdr=0.05,
                                        gp_version = GP_VERSION)
# write to out
fwrite(final_res_all,paste0("GO_res_Somatic_tissue_private_BND.txt"),row.names = F,na = NA,sep='\t',quote=F)



cat("\n################################# Running GO on all the somatic variants in simple repeats")
vcfbed_SR <- subset(vcfbed,final.anno=="simple_rep")
somatic_SR <- subset(vcfbed_SR,!is.na(vcfbed_SR$somatic.anno))
# run go
final_res_SR <- wrapGO_tissue_specific(background=somatic_SR,
                                       svtype="BND",
                                       svclass="Somatic_tissue_private",
                                       fdr=0.1,
                                       gp_version = GP_VERSION)
# write to out
fwrite(final_res_SR,paste0("GO_res_Simple_rep_Somatic_tissue_private_BND.txt"),row.names = F,na = NA,sep='\t',quote=F)



cat("\n################################# Running GO on all the somatic variants NOT in simple repeats")
vcfbed_nSR <- subset(vcfbed,final.anno!="simple_rep")
somatic_nSR <- subset(vcfbed_nSR,!is.na(vcfbed_nSR$somatic.anno))
# run go
final_res_nSR <- wrapGO_tissue_specific(background=somatic_nSR,
                                       svtype="BND",
                                       svclass="Somatic_tissue_private",
                                       fdr=0.1,
                                       gp_version = GP_VERSION)
# write to out
fwrite(final_res_nSR,paste0("GO_res_normal_Somatic_tissue_private_BND.txt"),row.names = F,na = NA,sep='\t',quote=F)

cat("Done\n")
cat("\n#################################################################################\n")
