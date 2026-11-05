###########################################################################################################
### 08.
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
GP_VERSION <- "e111_eg58_p18"   # results in this folder were produced with this g:Profiler release
gprofiler2::set_base_url(gp_archive_url(GP_VERSION))
vi <- tryCatch(gprofiler2::get_version_info("mmusculus"), error = function(e) NULL)
if (!is.null(vi)) cat("g:Profiler:", vi$gprofiler_version,
                      "|", vi$sources$`GO:BP`$version, "\n")
###########################################################################################################

cat("\n##################################################################################\n")
cat("Loading:\n\t../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz\n")

# Load vcf
vcfbed <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
# split annotations (BND have anno/anno annotation format)
vcfbed$annotation1 <- str_split_fixed(vcfbed$annotation,"/",2)[,1]
vcfbed$annotation2 <- str_split_fixed(vcfbed$annotation,"/",2)[,2]
vcfbed$gene_name1 <- str_split_fixed(vcfbed$gene_name,"/",2)[,1]
vcfbed$gene_name2 <- str_split_fixed(vcfbed$gene_name,"/",2)[,2]

# Load somatic obj
cat("Loading:\n\t../02_get_somatic/somatic.Rds\n")
obj <- readRDS("../02_get_somatic/somatic.Rds")
###########################################################################################################


###########################################################################################################

cat("\n################################# Running GO on all the somatic variants")
# background is any gene hosting a somatic variant
# (within the function only genes hosting genic variants are selected)
somatic_bg <- subset(vcfbed,!is.na(vcfbed$somatic.anno))



# 1. run go somatic common private
final_res_all_private <- wrapGO_common(background=somatic_bg,
                                       svtype="DEL",
                                       svclass="Somatic_common_private",
                                       fdr=0.05,
                                       gp_version=GP_VERSION)

# 2. run go somatic common recurrent
final_res_all_recurrent <- wrapGO_common(background=somatic_bg,
                                       svtype="DEL",
                                       svclass="Somatic_common_recurrent",
                                       fdr=0.05,
                                       gp_version=GP_VERSION)

# bind and write results
final_res_all_private$Type <- "Somatic_common_private"
final_res_all_recurrent$Type <- "Somatic_common_recurrent"
tot <- rbind(final_res_all_private,final_res_all_recurrent)
fwrite(tot,paste0("GO_res_Somatic_common_DEL.txt"),row.names = F,na = NA,sep='\t',quote=F)
###########################################################################################################

###########################################################################################################

cat("\n################################# Running GO on all the somatic variants in simple repeats")
vcfbed_SR <- subset(vcfbed,final.anno=="simple_rep")
somatic_SR <- subset(vcfbed_SR,!is.na(vcfbed_SR$somatic.anno))

# 3. run go somatic common private
final_res_all_private_SR <- wrapGO_common(background=somatic_SR,
                                          svtype="DEL",
                                          svclass="Somatic_common_private",
                                          fdr=0.1,
                                          gp_version=GP_VERSION)

# 4. run go somatic common recurrent
final_res_all_recurrent_SR <- wrapGO_common(background=somatic_SR,
                                          svtype="DEL",
                                          svclass="Somatic_common_recurrent",
                                          fdr=0.1,
                                          gp_version=GP_VERSION)

# bind and write results
final_res_all_private_SR$Type <- "Somatic_common_private"
final_res_all_recurrent_SR$Type <- "Somatic_common_recurrent"
tot_SR <- rbind(final_res_all_private_SR,final_res_all_recurrent_SR)
fwrite(tot_SR,paste0("GO_res_Simple_rep_Somatic_common_DEL.txt"),row.names = F,na = NA,sep='\t',quote=F)
###########################################################################################################

###########################################################################################################


cat("\n################################# Running GO on all the somatic variants NOT in simple repeats")
vcfbed_nSR <- subset(vcfbed,final.anno!="simple_rep")
somatic_nSR <- subset(vcfbed_nSR,!is.na(vcfbed_nSR$somatic.anno))


# 5. run go somatic common private
final_res_all_private_nSR <- wrapGO_common(background=somatic_nSR,
                                          svtype="DEL",
                                          svclass="Somatic_common_private",
                                          fdr=0.1,
                                          gp_version=GP_VERSION)

final_res_all_recurrent_nSR <- wrapGO_common(background=somatic_nSR,
                                            svtype="DEL",
                                            svclass="Somatic_common_recurrent",
                                            fdr=0.1,
                                            gp_version=GP_VERSION)
# final_res_all_recurrent_nSR is empty

# bind and write results
final_res_all_private_nSR$Type <- "Somatic_common_private"
tot_nSR <- final_res_all_private_nSR
fwrite(tot_nSR,paste0("GO_res_normal_Somatic_common_DEL.txt"),row.names = F,na = NA,sep='\t',quote=F)
###########################################################################################################

###########################################################################################################




cat("\n#################################################################################\n")
cat("Done\n")
cat("\n#################################################################################\n")



