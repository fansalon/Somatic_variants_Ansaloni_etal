###########################################################################################################
### 08.
### GO analyses on variants of interest - goseq
###########################################################################################################
rm(list = ls())
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



# get genes Somatic_common_private
qb_genes_pv <- getQueryGenes_common(background = subset(vcfbed,!is.na(vcfbed$somatic.anno)),
                                 svtype     = "DEL",
                                 svclass    = "Somatic_common_private")
# convert to lower, as getgo below returns capitalized gene names
qb_genes_pv$genes <- tolower(qb_genes_pv$genes) 


# get genes Somatic_common_private
qb_genes_rc <- getQueryGenes_common(background = subset(vcfbed,!is.na(vcfbed$somatic.anno)),
                                    svtype     = "DEL",
                                    svclass    = "Somatic_common_recurrent")
# convert to lower, as getgo below returns capitalized gene names
qb_genes_rc$genes <- tolower(qb_genes_rc$genes) 

qb_genes <- rbind(qb_genes_pv,qb_genes_rc)
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

# get background
bg_genes <- unique(subset(qb_genes,feat_go=="background")[,"genes"])

# gene->cat and cat->gene converter
go_ls    <- suppressMessages(getGOanno(bg_genes,"mm39","geneSymbol"))
gene2cat <- go_ls$gene_to_cat
cat2gene <- go_ls$cat_to_gene
###########################################################################################################

# get Type
types <- unique(qb_genes$somatic.anno)

# get comparisons
feat <- unique(qb_genes$feat_go)
feat <- feat[feat!="background"]


###########################################################################################################

# calc GO
tot_res <- data.frame()
for (tp in types) {
  cat(tp,"\n")
  for (ft in feat) {
    
    cat("\t",ft,"\n")
    
    # get genes of interest for this comparison
    d1 <- subset(qb_genes,somatic.anno==tp & feat_go==ft )
    
    # run analysis
    res <- runGOseq(d1,bg_genes,gene_l_anno_aggr,gene2cat,cat2gene,min_genes=round(5*nrow(d1)/100,0))

    # add labels
    res$somatic.anno <- tp
    res$features     <- ft
    
    # bind
    tot_res <- rbind(tot_res,res)
    
  }
}
###########################################################################################################

write.table(tot_res,"GO_res_Somatic_common_DEL.txt",quote=F,sep='\t',row.names = F)

cat("\n#################################################################################\n")
cat("Done\n")
cat("\n#################################################################################\n")



