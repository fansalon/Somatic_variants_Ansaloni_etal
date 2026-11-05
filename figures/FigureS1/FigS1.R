#######################################################################################################################################
### Figure S1 - sequencing coverage + total number of variants per sample
#######################################################################################################################################
rm(list = ls())
library(ggplot2)
library(data.table)
library(stringr)
library(ggpubr)
#######################################################################################################################################
df <- read.table("../data/mapping_stat.txt",header=T,sep='\t')
genome_l <- sum(read.table("../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr")[,2])
vcfbed <- fread("../../analyses/01_parse_vcf/parsed_vcf.bed.gz",sep='\t',header = T,data.table = F)
#######################################################################################################################################

### Panel A - sequencing coverage


# Calc coverage as mapped bases / total bases in the mm39 genome
df$coverage <- df$mapped_bases / genome_l

# Rep names
df$Rep <- paste0("Mouse ",df$Rep-1)

pA=ggplot(df,aes(x=Tissue,y=coverage,fill=Tissue)) +
  theme_classic() +
  geom_boxplot(alpha=.6,outliers = F) +
  geom_point(aes(shape=Rep),size=2,position = position_jitter()) +
  # x axis
  xlab("Tissues") +
  theme(axis.text.x=element_text(size = 10,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=12)) +
  # y axis
  ylab("Coverage (X)") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12))

#######################################################################################################################################

#######################################################################################################################################


### Panel B - raw number of variants per sample
ff <- data.frame()
for (sample in colnames(vcfbed[,grepl("_ALT",colnames(vcfbed))])) {
  if (!grepl("_THRE",sample)) {
    # count variants by SVTYPE
    d1 <- subset(vcfbed,vcfbed[,sample]>=1)
    count_sample <- as.data.frame(table(d1$SVTYPE))
    count_sample$sample <- sample
    # bind
    ff <- rbind(ff,count_sample)
  }
}

# modify sample name
ff$sample2 <- paste0(str_split_fixed(ff$sample,"_",4)[,3],"_",as.numeric(str_split_fixed(ff$sample,"_",4)[,2])-1)
ff$mouse <- paste0("Mouse ",str_split_fixed(ff$sample2,"_",2)[,2])
ff$tissue <- str_split_fixed(ff$sample2,"_",2)[,1]


pB=ggplot(ff,aes(x=tissue,y=Freq,fill=Var1)) +
  theme_classic() +
  geom_bar(stat="identity",color='black') +
  # annotation
 # geom_text(aes(label=Var1),position = position_stack(vjust = 0.5),size=2.7) +
  # x axis
  xlab("xx") +
  theme(axis.text.x=element_text(size = 9,angle = 30, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("Number of variants") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  # legend
  theme(legend.title = element_blank()) +
  # facet
  facet_grid(~mouse,scales = "free",space = "free") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold'))


#######################################################################################################################################

#######################################################################################################################################


pp=ggarrange(pA,pB,ncol=2,widths = c(.4,.6),labels = c("A","B"))

ggsave(filename='FigS1.pdf',plot=pp,width=230,height=80,units='mm')

