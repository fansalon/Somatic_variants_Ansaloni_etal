###########################################################################################################
### Parse germline variants from manta results - CTX_3 Illumina sequencing
###########################################################################################################
rm(list = ls())
suppressMessages(library(stringr))
suppressMessages(library(BiocParallel))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
###########################################################################################################
source("../../functions/functions_nanopore.R")
###########################################################################################################
inpvcf     <- "../../16_runManta_Illumina/output/germline/results/variants/diploidSV.vcf"
simple_rep <- "../../../zenodo/input_data/simpleRepeat.tot.bed"
tes        <- "../../../zenodo/input_data/rmsk.mm39.txt.filtered.bed"
###########################################################################################################


cat("###########################################################################################\n")
cat("Loading vcf file:",inpvcf,"\n")

# Load vcf
vcf <- read.table(inpvcf,comment.char = "#",sep='\t',header = F)
colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample")
                   
# Remove unused columns
vcf$REF <- NULL
###################################################################################################

### BND and other SVs are formatted differently
vcf.nonbnd <- vcf[!grepl("SVTYPE=BND",vcf$INFO),]
vcf.bnd <- vcf[grepl("SVTYPE=BND",vcf$INFO),]



######## Non BND
cat("\tFormatting non BND variants\n")

# SVLEN, SVTYPE and END
vcf.nonbnd$SVLEN <- abs(as.numeric(str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"SVLEN=",2)[,2],";",2)[,1]))
vcf.nonbnd$SVTYPE <- str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"SVTYPE=",2)[,2],";",2)[,1]
vcf.nonbnd$END <- as.numeric(str_split_fixed(str_split_fixed(vcf.nonbnd$INFO,"END=",2)[,2],";",2)[,1])
vcf.nonbnd$END <- ifelse(vcf.nonbnd$SVTYPE=="INS",vcf.nonbnd$POS+1,vcf.nonbnd$END)
vcf.nonbnd$breakpoint2 <- NA

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

# define if BND is inter or intrachromosomal
vcf.bnd$CHROM_2 <- str_split_fixed(vcf.bnd$breakpoint2,":",2)[,1]
vcf.bnd$SVTYPE <- ifelse(vcf.bnd$CHROM==vcf.bnd$CHROM_2,"BND_intrachrm",
                         ifelse(vcf.bnd$CHROM!=vcf.bnd$CHROM_2,"BND_interchrm","ERR"))

# Remove unused cols
vcf.bnd$ALT <- NULL
vcf.bnd$QUAL <- NULL
vcf.bnd$INFO <- NULL
vcf.bnd$CHROM_2 <- NULL
vcf.bnd$END <- NULL


# define if there are specular BND - ie bp1 of BNDx is equal to bp2 of BNDy
vcf.bnd$breakpoint1 <- paste0(vcf.bnd$CHROM,':',vcf.bnd$POS)
specular <- table(vcf.bnd$breakpoint1 %in% vcf.bnd$breakpoint2)['TRUE']
# if there are specular, remove them
if (specular>0) {
  cat(paste0("\t\tfound ",specular," / ",nrow(vcf.bnd)," specular BNDs\n"))
  # order
  for (i in 1:nrow(vcf.bnd)) {
    d1 <- vcf.bnd[i,]
    if (d1$breakpoint1 >= d1$breakpoint2) {
      chromosome  <- str_split_fixed(d1$breakpoint1,":",2)[,1]
      start       <- as.numeric(str_split_fixed(d1$breakpoint1,":",2)[,2])
      breakpoint2 <- d1$breakpoint2
    } else if (d1$breakpoint1 < d1$breakpoint2) {
      chromosome  <- str_split_fixed(d1$breakpoint2,":",2)[,1]
      start       <- as.numeric(str_split_fixed(d1$breakpoint2,":",2)[,2])
      breakpoint2 <- d1$breakpoint1
    } else {
      stop("Hey?!?")
    }
    # update values
    vcf.bnd[i,"CHROM"]       <- chromosome
    vcf.bnd[i,"POS"]         <- start
    vcf.bnd[i,"breakpoint2"] <- breakpoint2
  }
}

# now remove breakpoint1 and get unique list of BND
vcf.bnd$breakpoint1 <- NULL
vcf.bnd <- unique(vcf.bnd)

cat(paste0("\t\tnumber of BNDs reduced to: ",nrow(vcf.bnd),"\n"))


###################################################################################################
vcf.bnd$END <- vcf.bnd$POS + 1
vcf.bnd <- vcf.bnd[,colnames(vcf.nonbnd)]
vcf2 <- rbind(vcf.nonbnd,vcf.bnd)


# do a check
if (! table(is.na(vcf2[grepl("BND",vcf2$SVTYPE),"breakpoint2"]))["FALSE"]==nrow(vcf2[grepl("BND",vcf2$SVTYPE),])) {
  stop("ERROR: there are some BND not assigned to any breakpoint2")
}


# Get info on genotype and supporting reads
vcf2$GT <- str_split_fixed(vcf2$sample,":",6)[,1]

# if PR (supporting PE reads) are reported
vcf2$REF_PE_reads <- ifelse(grepl("PR",vcf2$FORMAT),
                            as.numeric(str_split_fixed(str_split_fixed(vcf2$sample,":",6)[,5],",",2)[,1]),
                            0)
vcf2$ALT_PE_reads <- ifelse(grepl("PR",vcf2$FORMAT),
                            as.numeric(str_split_fixed(str_split_fixed(vcf2$sample,":",6)[,5],",",2)[,2]),
                            0)

# same for SR (supporting split reads)
vcf2$REF_split_reads <- ifelse(grepl("SR",vcf2$FORMAT),
                               as.numeric(str_split_fixed(str_split_fixed(vcf2$sample,":",6)[,6],",",2)[,1]),
                               0)
vcf2$ALT_split_reads <- ifelse(grepl("SR",vcf2$FORMAT),
                               as.numeric(str_split_fixed(str_split_fixed(vcf2$sample,":",6)[,6],",",2)[,2]),
                               0)
vcf2$REF_reads <- vcf2$REF_PE_reads + vcf2$REF_split_reads
vcf2$ALT_reads <- vcf2$ALT_PE_reads + vcf2$ALT_split_reads
vcf2$sample <- NULL
vcf2$MAF <- vcf2$ALT_reads / (vcf2$ALT_reads+vcf2$REF_reads)

# Convert to bed
vcf2$SV_ID <- paste(vcf2$CHROM,vcf2$POS,vcf2$END,vcf2$SVTYPE,vcf2$SVLEN,sep=":")
bed <- vcf2[,c("CHROM","POS","END","ID","SV_ID","FILTER","SVLEN","SVTYPE","breakpoint2","GT",
               "REF_PE_reads","ALT_PE_reads","REF_split_reads","ALT_split_reads","REF_reads","ALT_reads","MAF")]


###################################################################################################


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
# # Select the lucky ones
# dupl_selected <- as.data.frame(dupl %>% 
#                                  group_by(SV_ID) %>% 
#                                  slice_max(ALT_reads, n = 1,with_ties=F) )
# # Remove duplicated from original file
# n1 <- length(unique(bed$SV_ID))
# bed <- subset(bed,is.dupl==FALSE)
# # Add the selected duplicated
# bed <- rbind(bed,dupl_selected)
# rownames(bed) <- 1:nrow(bed)
# # Sort
# bed <- bed[order(bed$CHROM,bed$POS),]
# # Check
# n2 <- length(unique(bed$SV_ID))
# if (n1!=n2) { stop("Some SVs have been lost") }
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


###################################################################################################


#### Overlap with simple rep
cat("Overlapping with Simple repeats\n")
write.table(bed,"bedA.bed",row.names = F,col.names = F,sep='\t',quote=F)
system(paste0("bedtools intersect -nonamecheck -c -f .5 -a bedA.bed -b ",simple_rep," > intersect"))

# read file
res_int <- fread("intersect",data.table = F)
system("rm intersect bedA.bed")

# get SV ID overlapping simple repeats
res_int_rep <- subset(res_int,V19>=1)
black_list <- unique(res_int_rep$V5)

# classify
if (table(black_list %in% bed$SV_ID)["TRUE"]!=length(black_list)) {stop("Some of the SVs in the black list are missing from the tot file")}
bed$overlaps.simplerep <- bed$SV_ID %in% black_list
cat("\t",length(black_list),"variants overlapping simple repeats\n")


##########################################################################################


#### Overlap with TEs
cat("Overlapping with TEs\n")
write.table(bed,"bedA.bed",row.names = F,col.names = F,sep='\t',quote=F)
system(paste0("bedtools intersect -nonamecheck -c -f .5 -a bedA.bed -b ",tes," > intersect"))

# read file
res_int <- fread("intersect",data.table = F)
system("rm intersect bedA.bed")

# get SV ID overlapping simple repeats
res_int_rep <- subset(res_int,V20>=1)
black_list <- unique(res_int_rep$V5)

# classify
if (table(black_list %in% bed$SV_ID)["TRUE"]!=length(black_list)) {stop("Some of the SVs in the black list are missing from the tot file")}
bed$overlaps.TEs <- bed$SV_ID %in% black_list

cat("\t",length(black_list),"variants overlapping TEs\n")

###################################################################################################

# rename ID for BND
bed$SV_ID <- ifelse(grepl("BND",bed$SVTYPE),paste0(bed$CHROM,":",bed$POS,"@",bed$breakpoint2,":",bed$SVTYPE),bed$SV_ID)


##########################################################################################

### write to out

colnames(bed)[1] <- paste0("#",colnames(bed)[1])

write.table(bed,paste0(inpvcf,".bed"),row.names = F, quote = F, sep = '\t')
system(paste0("mv ",inpvcf,".bed ../"))
cat("###########################################################################################\n")
