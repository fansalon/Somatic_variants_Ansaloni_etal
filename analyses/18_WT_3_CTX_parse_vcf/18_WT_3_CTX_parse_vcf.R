###########################################################################################################
### Parse WT_3_CTX variants from sniffels
###########################################################################################################
rm(list = ls())
suppressMessages(library(stringr))
suppressMessages(library(BiocParallel))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
###########################################################################################################
source("../functions/functions_nanopore.R")
###########################################################################################################
inpvcf     <- "../../zenodo/input_data/WT_3_CTX.vcf"
simple_rep <- "../../zenodo/input_data/simpleRepeat.tot.bed"
tes        <- "../../zenodo/input_data/rmsk.mm39.txt.filtered.bed"
###########################################################################################################


cat("###########################################################################################\n")
cat("Loading vcf file:",inpvcf,"\n")

# Load vcf
vcf <- read.table(inpvcf,comment.char = "#",sep='\t',header = F)
colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample")
                   
# Remove unused columns
vcf$REF <- NULL
vcf$QUAL <- NULL
vcf$FILTER <- NULL
###################################################################################################

### BND and other SVs are formatted differently
vcf.nonbnd <- vcf[!grepl("SVTYPE=BND",vcf$INFO),]
vcf.bnd <- vcf[grepl("SVTYPE=BND",vcf$INFO),]



######## Non BND
cat("\tFormatting non BND variants\n")

# SVLEN, SVTYPE, END, 2nd breakp., STDEV_LEN and STDEV_POS, PRECISE, MOSAIC
vcf.nonbnd$SVLEN <- abs(as.numeric(str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"SVLEN=",2)[,2],";",2)[,1]))
vcf.nonbnd$SVTYPE <- str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"SVTYPE=",2)[,2],";",2)[,1]
vcf.nonbnd$END <- as.numeric(str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"END=",2)[,2],";",2)[,1])
vcf.nonbnd$END <- ifelse(vcf.nonbnd$SVTYPE=="INS",vcf.nonbnd$POS+1,vcf.nonbnd$END)
vcf.nonbnd$breakpoint2 <- NA
vcf.nonbnd$STDEV_LEN <- as.numeric(str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"STDEV_LEN=",2)[,2],";",2)[,1])
vcf.nonbnd$STDEV_POS <- as.numeric(str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"STDEV_POS=",2)[,2],";",2)[,1])
vcf.nonbnd$PRECISE <- str_split_fixed(vcf.nonbnd$INFO,";",2)[,1]
vcf.nonbnd$MOSAIC <- ifelse(grepl("MOSAIC",vcf.nonbnd$INFO),"MOSAIC","nonMOSAIC")


# Remove unused cols
vcf.nonbnd$ALT <- NULL
vcf.nonbnd$QUAL <- NULL
vcf.nonbnd$INFO <- NULL


######## BND

cat("\tFormatting BND variants\n")

# SVLEN, SVTYPE and END
vcf.bnd$SVLEN <- 1
vcf.bnd$SVTYPE <- str_split_fixed(str_split_fixed(vcf.bnd$INFO,"SVTYPE=",2)[,2],";",2)[,1]
vcf.bnd$END <- vcf.bnd$POS + 1
# Get second breakpoint for BND
vcf.bnd$breakpoint2 <- vcf.bnd$ALT
# remove "[" or "]" from breakpoint
vcf.bnd$breakpoint2 <- gsub("\\[","@",vcf.bnd$breakpoint2)
vcf.bnd$breakpoint2 <- gsub("\\]","@",vcf.bnd$breakpoint2)
vcf.bnd$breakpoint2 <- str_split_fixed(vcf.bnd$breakpoint2,"@",3)[,2]


# get STDEV_LEN and STDEV_POS, PRECISE, MOSAIC
vcf.bnd$STDEV_LEN <- 0
vcf.bnd$STDEV_POS <- as.numeric(str_split_fixed(str_split_fixed(vcf.bnd$INFO,"STDEV_POS=",2)[,2],";",2)[,1])
vcf.bnd$PRECISE <- str_split_fixed(vcf.bnd$INFO,";",2)[,1]
vcf.bnd$MOSAIC <- ifelse(grepl("MOSAIC",vcf.bnd$INFO),"MOSAIC","nonMOSAIC")


# Remove unused cols
vcf.bnd$ALT <- NULL
vcf.bnd$QUAL <- NULL
vcf.bnd$INFO <- NULL
vcf.bnd$CHROM_2 <- NULL

###################################################################################################

vcf2 <- rbind(vcf.nonbnd,vcf.bnd)

# do a check
if (! table(is.na(vcf2[grepl("BND",vcf2$SVTYPE),"breakpoint2"]))["FALSE"]==nrow(vcf2[grepl("BND",vcf2$SVTYPE),])) {
  stop("ERROR: there are some BND not assigned to any breakpoint2")
}


# Get info on genotype and supporting reads
vcf2$GT <- str_split_fixed(vcf2$sample,":",4)[,1]


vcf2$REF_reads <- as.numeric(str_split_fixed(vcf2$sample,":",4)[,3])
vcf2$ALT_reads <- as.numeric(str_split_fixed(vcf2$sample,":",4)[,4])
vcf2$sample <- NULL
vcf2$FORMAT <- NULL
vcf2$MAF <- vcf2$ALT_reads / (vcf2$ALT_reads+vcf2$REF_reads)

# Convert to bed
vcf2$SV_ID <- paste(vcf2$CHROM,vcf2$POS,vcf2$END,vcf2$SVTYPE,vcf2$SVLEN,sep=":")
bed <- vcf2[,c("CHROM","POS","END","ID","SV_ID","SVLEN","STDEV_LEN","STDEV_POS","SVTYPE","breakpoint2",
               "PRECISE","MOSAIC","GT","REF_reads","ALT_reads","MAF")]


###################################################################################################
# 
# 
# ### Remove duplicated SVs
# 
# cat(paste0("\tSorting duplicated variants\n\t\tstarting from         : ",nrow(bed),"\n"))
# 
# # for BND generate an ID that takes into account both breakpoints, 
# # ie, two BND can have the same bp#1, but different bp#2 
# # in these cases, I want them to be considered as independent events
# bed$SV_ID <- ifelse(grepl("BND",bed$SVTYPE),paste(bed$CHROM,bed$POS,gsub(":","_",bed$breakpoint2),bed$SVTYPE,"1",sep=":"),bed$SV_ID)
# 
# ### Among the duplicated SVs, select the ones with the highest mean of alternative reads
# # Mark dupl
# dupl.id <- bed[ duplicated(bed$SV_ID),"SV_ID"]
# bed$is.dupl <- bed$SV_ID %in% dupl.id
# # Get dupl. SVs
# dupl <- bed[bed$is.dupl,]
# 
# if (nrow(dupl)>0) {
# 
#   # Select the lucky ones
#   dupl_selected <- as.data.frame(dupl %>% 
#                                    group_by(SV_ID) %>% 
#                                    slice_max(ALT_reads, n = 1,with_ties=F) )
#   # Remove duplicated from original file
#   n1 <- length(unique(bed$SV_ID))
#   bed <- subset(bed,is.dupl==FALSE)
#   # Add the selected duplicated
#   bed <- rbind(bed,dupl_selected)
#   rownames(bed) <- 1:nrow(bed)
#   # Sort
#   bed <- bed[order(bed$CHROM,bed$POS),]
#   # Check
#   n2 <- length(unique(bed$SV_ID))
#   if (n1!=n2) { stop("Some SVs have been lost") }
# }
# 
# 
# cat(paste0("\t\tduplicated variants   : ",nrow(dupl),"\n"))
# cat(paste0("\t\tnon-redundant variants: ",nrow(bed),"\n"))

bed$is.dupl <- NA


###################################################################################################

# cat(paste0("\tRemoving variants on scaffolds or on chrM\n"))
# 
# 
# nn1 <- nrow(bed)
# 
# ### Remove SV on scaffolds
# bed <- bed[grepl("chr",bed$CHROM),]
# # discard scaffolds also from breakpoint2 in BND
# bed <- bed[ is.na(bed$breakpoint2) | grepl("chr",bed$breakpoint2),]
# # remove also chrM
# bed <- bed[!grepl("chrM",bed$breakpoint2),]
# 
# nn2 <- nrow(bed)
# 
# cat(paste0("\t\tremoved ",nn1-nn2," variants\n"))



##########################################################################################


cat("Overlapping with Simple repeats\n")


write.table(bed,"bedA.bed",row.names = F,col.names = F,sep='\t',quote=F)
system(paste0("bedtools intersect -nonamecheck -c -f .5 -a bedA.bed -b ",simple_rep," > intersect"))

# read file
res_int <- fread("intersect",data.table = F)
system("rm intersect bedA.bed")

# get SV ID overlapping simple repeats
res_int_rep <- subset(res_int,V18>=1)
black_list <- unique(res_int_rep$V5)

# classify
if (table(black_list %in% bed$SV_ID)["TRUE"]!=length(black_list)) {stop("Some of the SVs in the black list are missing from the bed file")}
bed$overlaps.simplerep <- bed$SV_ID %in% black_list
cat("\t",length(black_list),"variants overlapping simple repeats\n")

###################################################################################################


cat("Overlapping with TEs\n")


write.table(bed,"bedA.bed",row.names = F,col.names = F,sep='\t',quote=F)
system(paste0("bedtools intersect -nonamecheck -c -f .5 -a bedA.bed -b ",tes," > intersect"))

# read file
res_int <- fread("intersect",data.table = F)
system("rm intersect bedA.bed")

# get SV ID overlapping simple repeats
res_int_rep <- subset(res_int,V19>=1)
black_list <- unique(res_int_rep$V5)

# classify
if (table(black_list %in% bed$SV_ID)["TRUE"]!=length(black_list)) {stop("Some of the SVs in the black list are missing from the bed file")}
bed$overlaps.TEs <- bed$SV_ID %in% black_list

cat("\t",length(black_list),"variants overlapping TEs\n")

###################################################################################################

# rename ID for BND
bed$SV_ID <- ifelse(grepl("BND",bed$SVTYPE),paste0(bed$CHROM,":",bed$POS,"@",bed$breakpoint2,":",bed$SVTYPE),bed$SV_ID)


###################################################################################################




###################################################################################################

### write to out

colnames(bed)[1] <- paste0("#",colnames(bed)[1])

write.table(bed,paste0(inpvcf,".bed"),row.names = F, quote = F, sep = '\t')
system(paste0("mv ",inpvcf,".bed ."))
cat("###########################################################################################\n")
