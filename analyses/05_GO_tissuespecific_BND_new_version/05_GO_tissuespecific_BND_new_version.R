###########################################################################################################
### 05.
### This script performs GO analysis on BND
###########################################################################################################
rm(list = ls())
seed_f <- "11322"
set.seed(seed_f)
###########################################################################################################
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(gprofiler2))
suppressMessages(library(rtracklayer) )
suppressMessages(library(goseq))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################
anno <- "../../zenodo/input_data/gencode.GRCm39.vM28.annotation.gtf"
###########################################################################################################

cat("\n##################################################################################\n")

# Load gtf
cat("Loading:\n\t../../zenodo/input_data/gencode.GRCm39.vM28.annotation.gtf\n")
gtf <- import(paste0(anno))

# Load vcf
cat("\t../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz\n")
vcfbed <- fread("../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
# split annotations
vcfbed$annotation1 <- str_split_fixed(vcfbed$annotation,"/",2)[,1]
vcfbed$annotation2 <- str_split_fixed(vcfbed$annotation,"/",2)[,2]
vcfbed$gene_name1 <- str_split_fixed(vcfbed$gene_name,"/",2)[,1]
vcfbed$gene_name2 <- str_split_fixed(vcfbed$gene_name,"/",2)[,2]

# Load somatic obj
cat("\t../02_get_somatic/somatic.Rds\n")
obj <- readRDS("../02_get_somatic/somatic.Rds")
###########################################################################################################



###########################################################################################################
# get genes Somatic_tissue_private/BND
qb_genes <- getQueryGenes(background = subset(vcfbed,!is.na(vcfbed$somatic.anno)),
                         svtype     = "BND",
                         svclass    = "Somatic_tissue_private")
# convert to lower, as getgo below returns capitalized gene names
qb_genes$genes <- tolower(qb_genes$genes) 
###########################################################################################################



###########################################################################################################
# get gene length
gtf_g            <- as.data.frame(subset(gtf,type=="gene"))
gene_l_anno_aggr <- getGeneLength(gtf_g)
###########################################################################################################



###########################################################################################################
# Since I have converted all the genes to lower case, I keep track of the conversion to the original
# gene name so I can then convert this back
sym_map <- setNames(gtf_g$gene_name, tolower(gtf_g$gene_name))
sym_map <- sym_map[!duplicated(names(sym_map))]
###########################################################################################################



###########################################################################################################
# get background genes and gene->cat and cat->gene converter

# get background
bg_genes <- unique(subset(qb_genes,Type=="background")[,"genes"])

go_ls    <- suppressMessages(getGOanno(bg_genes,"mm39","geneSymbol"))
gene2cat <- go_ls$gene_to_cat
cat2gene <- go_ls$cat_to_gene
###########################################################################################################



# get Type
types <- unique(qb_genes$Type)
types <- types[types!="background"]

# get comparisons
feat <- unique(qb_genes$feat_go)
feat <- feat[feat!="background"]



# calc GO
tot_res <- data.frame()
for (tp in types) {
  cat(tp,"\n")
  for (ft in feat) {
    
    cat("\t",ft,"\n")
    
    # get genes of interest for this comparison
    d1 <- subset(qb_genes,Type==tp & feat_go==ft )
    
    # run analysis
    res <- runGOseq(d1,bg_genes,gene_l_anno_aggr,gene2cat,cat2gene,min_genes=round(5*nrow(d1)/100,0))
    
    
    # add labels
    res$Type     <- tp
    res$features <- ft
    
    # bind
    tot_res <- rbind(tot_res,res)
    
  }
}

write.table(tot_res,"GO_res_Somatic_tissue_private_BND.txt",quote=F,sep='\t',row.names = F)




