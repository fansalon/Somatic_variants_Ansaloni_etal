library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(stringr)

bed="../../03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz"
sr1="../../../zenodo/input_data/simpleRepeat.rmsk.bed"
sr2="../../../zenodo/input_data/simpleRepeat.txt.bed"


# 1 get coords variants entirely contained in SR
df=fread(bed,header=T,sep='\t',data.table=F)
df=subset(df,SVTYPE=="DEL")
df2=subset(df,intersect_SR_perc>=100)[,c("chrom","start","end","SVTYPE","SVLEN","somatic.anno","Samples","SR","SR_ID")]
cat(paste0("Selected ",nrow(df2)," DELs entirely contained in simple repeats\n"))
# get somatic only
df2=subset(df2,!is.na(df2$somatic.anno))
cat(paste0("\tselected ",nrow(df2)," somatic DELs only\n"))
# get somatic of interest
df2=subset(df2,somatic.anno=="Somatic_common_private" |
             somatic.anno=="Somatic_common_recurrent")
cat(paste0("\tselected ",nrow(df2)," somatic DELs Somatic_common_private or Somatic_common_recurrent\n"))
(table(df2$somatic.anno))
nro1=nrow(df2)


# 2 get coords of the SR themselves
df2$id <- paste(df2$chrom,df2$start,df2$end,df2$SVTYPE,df2$SVLEN,df2$somatic.anno,df2$Samples,sep = "-")
df2 <- df2[,c("SR_ID","id","SR")]
df2$chrom <- str_split_fixed(df2$SR_ID,"-",3)[,1]
df2$start <- as.numeric(str_split_fixed(df2$SR_ID,"-",3)[,2])
df2$end <- as.numeric(str_split_fixed(df2$SR_ID,"-",3)[,3])
df3 <- df2[,c("chrom","start","end","SR","id")]


# 3 format as requested by nanorepeat
#delet (, ) and n
df3$SR <- gsub("\\(","",df3$SR)
df3$SR <- gsub("\\)n","",df3$SR)


# 5 look at motif and SR length distribution
df3$motif_l <- nchar(df3$SR)
df3$len <- df3$end - df3$start
len_95 <- quantile(df3$len,probs = .95)
len_motif_95 <- quantile(df3$motif_l,probs = .95)

p1=ggplot(df3) +
  theme_classic() +
  geom_histogram(aes(len)) +
  # x axis
  xlab("simple repeat length [bp]") +
  theme(axis.text.x=element_text(size = 10,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=12)) +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  # y axis
  ylab("No. of SR hosting\nsomatic common DELs") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12))  +
  annotate("text",x=2500,y=2000,label=paste0("95th percentile = ",round(len_95,0)),
           size=3,color='red') +
  geom_vline(xintercept = len_95,
             linetype='dashed',color='red')
p2=ggplot(df3) +
  theme_classic() +
  geom_histogram(aes(motif_l)) +
  # x axis
  xlab("simple repeat motif length [bp]") +
  theme(axis.text.x=element_text(size = 10,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=12)) +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  # y axis
  ylab("No. of SR hosting\nsomatic common DELs") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  annotate("text",x=85,y=5400,label=paste0("95th percentile = ",len_motif_95),
           size=3,color='red') +
  geom_vline(xintercept = len_motif_95,
             linetype='dashed',color='red')

ff=ggarrange(p1,p2,ncol=2)
#ggsave(filename='plot.pdf',plot=ff,width=210,height=70,units='mm')


# discard SR longer than 95% and with motif longer than 95%
df3 <- subset(df3, len < len_95 &
                     motif_l < len_motif_95 )
df3$motif_l=NULL
df3$len=NULL
cat(paste0("\tselected ",nrow(df3),
           " simple repeats hosting somatic DELs Somatic_common_private or Somatic_common_recurrent shorter than ",
           round(len_95,0),"bp (95th percentile) and with a SR motif shorter than ",len_motif_95,"bp (95th percentile)\n"))



# write to out
fwrite(df3,"input_SR.bed",col.names=F,row.names=F,sep='\t',na=NA,quote=F)
