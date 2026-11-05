###########################################################################################################
### 15.
### GO analyses on MEI
###########################################################################################################
rm(list = ls())
###########################################################################################################
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(gprofiler2))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################

cat("\n##################################################################################\n")

# Load vcf - for background
vcfbed <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
mei <- read.table("../14_annotate_MEI/genomic_feature_anno.txt",header=T,sep='\t')
###########################################################################################################
fdr=0.9
fdr=1

###########################################################################################################

cat("\n################################# Running GO on all the somatic MEI")


# 1 get background
# background is any gene hosting a somatic variant
cat("\nGetting background genes (all the genes that host a genic somatic variant)")
background <- subset(vcfbed,!is.na(vcfbed$somatic.anno))

# split annotations (BND have anno/anno annotation format)
background$annotation1 <- str_split_fixed(background$annotation,"/",2)[,1]
background$annotation2 <- str_split_fixed(background$annotation,"/",2)[,2]
background$gene_name1 <- str_split_fixed(background$gene_name,"/",2)[,1]
background$gene_name2 <- str_split_fixed(background$gene_name,"/",2)[,2]

# Gene name relative to genes hosting somatic variants
all_genes_genic_variant1 <- unique(subset(background,annotation1!="Distal Intergenic" &
                                            annotation1!="Downstream")[,"gene_name1"])
all_genes_genic_variant2 <- unique(subset(background,annotation2!="Distal Intergenic" &
                                            annotation2!="Downstream")[,"gene_name2"])
all_genes_genic_variant <- unique(c(all_genes_genic_variant1,all_genes_genic_variant2))
# discard empty (not all the variants have a 2nd gene, but only the BND)
all_genes_genic_variant <- all_genes_genic_variant[!all_genes_genic_variant==""]
cat(": ",length(all_genes_genic_variant),"\n")



# 2 get somatic MEI
cat(paste0("Getting somatic MEI (",nrow(mei),")\n"))


### START ANALYSIS


###### Start GO

# Get columns of interest (annotation and gene name)
input_df <- unique(mei[,c("ID","annotation","gene_name")])


# Comparison by comparison GO analysis
features <- c("Promoter","Intron")
counter <- 0
for (feat in features) {
  cat("\t\t",feat,":")
  
  # Get genes associated to the feature
  d1_genes <- unique(subset(input_df,annotation==feat)[,"gene_name"])
  cat(length(d1_genes),"\n")
  # run go
  n_significant_genes <- round(5*length(d1_genes)/100,0)
  cat("\t\tGO must be found in at least 5% of the genes (",n_significant_genes,"genes)\n")
  cat("\t\tFDR =",fdr,"\n")
  res <- doGprof2(query_gene = d1_genes,bg_genes = all_genes_genic_variant,
                  fdr = fdr,GO_name = "GO:BP",n_significant_genes = n_significant_genes)
  # bind
  if (nrow(res)>0) {
    counter <- counter + 1
    res$annotation <- feat
    if (counter==1) {
      final_res <- res
    } else {
      final_res <- rbind(final_res,res)
    }
  }
}
if (!exists("final_res")) {final_res=data.frame()}
fwrite(final_res,paste0("GO_res.txt"),row.names = F,na = NA,sep='\t',quote=F)


cat("\n#################################################################################\n")
cat("Done\n")
cat("\n#################################################################################\n")



