#######################################################################################################################################
### Figure 1
#######################################################################################################################################
rm(list = ls())
library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(ggpubr)
library(circlize)
#######################################################################################################################################
### Graphical parameters
# ggplot
alpha=.8
axistexx=10
axistitlex=12
axistexty=10
axistitley=12
#######################################################################################################################################
#setwd("~/Dropbox/IIT/WT_Mouse/paper/figures/Figure1/")
vcf.final <- fread("../../analyses/01_parse_vcf/parsed_vcf.bed.gz",header=T,sep='\t',data.table = F)
rownames(vcf.final) <- vcf.final$ID
obj <- readRDS("../../analyses/02_get_somatic/somatic.Rds")
#######################################################################################################################################


#######################################################################################################################################
####### Panel B - circos

# read bed file with chrom length and format it
cytoband.file="../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr"
cytoband.df = read.table(cytoband.file,  sep = "\t")
# convert it to bed
names(cytoband.df) <- c("chr","end")
cytoband.df$start <- 1
cytoband.df <- cytoband.df[,c("chr","start","end")]
# discard chrM
cytoband.df <- subset(cytoband.df,chr!='chrM')
cytoband.df$V4 <- "."
cytoband.df$V5 <- "."

# total number of variants, reported with , as thousands separator
name=paste0(gsub("\\.",",",nrow(vcf.final)/1000),"\ncombined variants")
colo="orange"

# Select SV that are not BND
no_bnd <- subset(vcf.final,SVTYPE!="BND")[,c("chrom","start","end")]
names(no_bnd)[1] <- "chr"
dels <- no_bnd[grepl("DEL",rownames(no_bnd)),]
dup <- no_bnd[grepl("DUP",rownames(no_bnd)),]
ins <- no_bnd[grepl("INS",rownames(no_bnd)),]

# Select BND only and get first breakpoint only
bnd <- subset(vcf.final,SVTYPE=="BND")[,c("chrom","start","end","breakpoint2")]
names(bnd)[1] <- "chr"
bnd1 <- bnd[,c("chr","start","end")]


pdf('PanelB.pdf', width=4.5, height=4.5)
circos.clear()
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.05)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
circos.initializeWithIdeogram(cytoband.df,plotType=c("labels"))
#circos.genomicDensity(inp_ls, col = colo,count_by="number",track.height=.2)
circos.genomicDensity((dels), col = "orange",count_by="number",track.height=.1,window.size = 3e6)
circos.genomicDensity((ins), col = "yellow3",count_by="number",track.height=.1,window.size = 8e6)
circos.genomicDensity((dup), col = "lightblue",count_by="number",track.height=.1,window.size = 6e6)
circos.genomicDensity((bnd1), col = "dodgerblue",count_by="number",track.height=.1,window.size = 5e6)
#circos.genomicLink(head(bnd1,10000), head(bnd2,10000),col="grey90",lwd=0.0000000000000000000000001)
text(0, 0, name, cex = .79, font=2)
# Add a legend
legend(.9, -.8, legend=c("DEL","INS","DUP","BND"),col=c("orange","yellow3","lightblue","dodgerblue"), pch=16, cex=0.8)
dev.off()


#######################################################################################################################################

######### Panel C: classify the variants as 0/0, 0/1 and 1/1
genotype_count <- obj$GT_count

panlC=ggplot(genotype_count,aes(x="Combined\nvariants",y=Freq,fill=Var1)) +
  theme_classic() +
  geom_bar(stat="identity",color='black',aes(fill=(Var1)),alpha=alpha) +
  # annotation
  geom_text(aes(label=Var1),position = position_stack(vjust = 0.5),size=2.7) +
  #annotate("text",x=1.9,y=160000,label=paste0("Combined variants: ",tot)) +
  # x axis
  xlab("Supporting reads") +
  theme(axis.text.x=element_text(size = axistexx,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("No. of combined variants") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=axistitley)) +
  # legend
  theme(legend.position ="none") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("0/0"="seagreen3",
                               "0/1"="deepskyblue1",
                               "1/1"="lightpink1")) +
  annotate("text",x=1,y=50000,label=subset(genotype_count,Var1=="0/0")[,"Freq"],size=2.7)

#######################################################################################################################################


#######################################################################################################################################

### Panel D
final_count.aggr <- obj$somaticSV_mouse_tissue_count

panelD=ggplot(final_count.aggr,aes(x=Var1,y=Freq)) +
  theme_classic() +
  geom_bar(stat="identity",aes(fill=Var1),color='black',alpha=alpha) +
  # facet
  facet_grid(~Mouse,scales = "free",space = "free") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # annotation
  geom_text(aes(y=Freq+1500,label=paste0(round(Freq/Tot*100,0),"%")),size=3) +
  geom_text(aes(y=32000,label=paste0("n = ",Tot),x=1.8),size=2.8) +
  # x axis
  xlab("Number of tissues sharing the 0/0 variants") +
  theme(axis.text.x=element_text(size = axistexx,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=axistitlex)) +
  # y axis
  ylab("No. of 0/0 variants") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=axistitley)) +
  # legend
  theme(legend.position ="none") +
  theme(legend.title = element_blank()) +
  # scale fill
  scale_fill_manual(values = c("1/3"=hue_pal()(3)[2],"2/3"=hue_pal()(3)[2],"3/3"="seagreen4"))
#######################################################################################################################################



#######################################################################################################################################


#### Panel E - circos 3/3 variants
common_final <- obj$common_somatic
# 33621

### Draw circos plot with these common somatic variants

# Get coords
common_final_bed <- data.frame(chr=str_split_fixed(rownames(common_final),"-",4)[,1],
                               start=as.numeric(str_split_fixed(rownames(common_final),"-",4)[,2]),
                               end=as.numeric(str_split_fixed(rownames(common_final),"-",4)[,3]),
                               svtype=str_split_fixed(rownames(common_final),"-",4)[,4])

# Select SV that are not BND
common_final_DEL <- subset(common_final_bed,svtype=="DEL")
common_final_DUP <- subset(common_final_bed,svtype=="DUP")
common_final_INS <- subset(common_final_bed,svtype=="INS")

# Select BND only and get first breakpoint only
common_final_BND <- subset(common_final_bed,svtype=="BND")

name=paste0(gsub("\\.",",",nrow(common_final)/1000),"\ncommon\n(3/3 tissues)\nsomatic variants")
pdf('PanelE.pdf', width=4.5, height=4.5)
circos.clear()
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.05)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
circos.initializeWithIdeogram(cytoband.df,plotType=c("labels"))
circos.genomicDensity(common_final_DEL, col = "orange",count_by="number",track.height=.1,
                      window.size = 3e6)
circos.genomicDensity(common_final_INS, col = "yellow3",count_by="number",track.height=.1,
                      window.size = 3e6)
circos.genomicDensity(common_final_DUP, col = "lightblue",count_by="number",track.height=.1,
                      window.size = 3e6)
circos.genomicDensity(common_final_BND, col = "dodgerblue",count_by="number",track.height=.1,
                      window.size = 3e6)
text(0, 0, name, cex = .79, font=2)
dev.off()
#######################################################################################################################################


#######################################################################################################################################

#### Panel F - circos 1/3 and 2/3 variants
private_final <- obj$tissue_specific_somatic
nrow(private_final)
# 82641


### Draw circos plot with these common somatic variants

# Get coords
private_final_bed <- data.frame(chr=str_split_fixed(rownames(private_final),"-",4)[,1],
                                start=as.numeric(str_split_fixed(rownames(private_final),"-",4)[,2]),
                                end=as.numeric(str_split_fixed(rownames(private_final),"-",4)[,3]),
                                svtype=str_split_fixed(rownames(private_final),"-",4)[,4])

# Select SV that are not BND
private_final_DEL <- subset(private_final_bed,svtype=="DEL")
private_final_DUP <- subset(private_final_bed,svtype=="DUP")
private_final_INS <- subset(private_final_bed,svtype=="INS")

# Select BND only and get first breakpoint only
private_final_BND <- subset(private_final_bed,svtype=="BND")


name=paste0(gsub("\\.",",",nrow(private_final_bed)/1000),"\ntissue-specific\n(1/3 and 2/3 tissues)\nsomatic variants")
pdf('PanelF.pdf', width=4.5, height=4.5)
circos.clear()
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.05)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
circos.initializeWithIdeogram(cytoband.df,plotType=c("labels"))
circos.genomicDensity(private_final_DEL, col = "orange",count_by="number",track.height=.1,
                      window.size = 3e6)
circos.genomicDensity(private_final_INS, col = "yellow3",count_by="number",track.height=.1,
                      window.size = 3e6)
circos.genomicDensity(private_final_DUP, col = "lightblue",count_by="number",track.height=.1,
                      window.size = 4e7)
circos.genomicDensity(private_final_BND, col = "dodgerblue",count_by="number",track.height=.1,
                      window.size = 3e6)
text(0, 0, name, cex = .79, font=2)
dev.off()
#######################################################################################################################################


riga2=ggarrange(NULL,panlC,panelD,widths = c(.38,.17,.55),labels = c("BA","C","DEF"),nrow=1)
ff=ggarrange(NULL,riga2,nrow=3)

ggsave(filename='PanelCD.pdf',plot=ff,width=230,height=250/4*3.5,units='mm')


