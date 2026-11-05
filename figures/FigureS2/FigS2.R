#######################################################################################################################################
### Figure S2 - distribution genomic features - tissue specific
###############################################################################################
library(ggplot2)
library(reshape2)
library(ggpubr)
###############################################################################################
downsample   <- read.table("../../analyses/22_downsampling/output/downsampling_WT_3_HIP/1/variant_count.txt",header=T)
tot_variants <- read.table("../../analyses/22_downsampling/output/count_variants_allsamples/variant_count_ALL.txt",header=T)
genome_l     <- sum(read.table("../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr")[,2])
mapped       <- read.table("../data/mapping_stat.txt",header=T) 
###############################################################################################



###############################################################################################
### Panel A - downsample of HIP-3

pvalue_all <- formatC(cor.test(downsample$Variants,downsample$Mapped_reads)$p.value,digits=2,format='e')
pvalue_smt <- formatC(cor.test(downsample$Somatic,downsample$Mapped_reads)$p.value,digits=2,format='e')

df2 <- downsample
df2$Downsample <- paste0(df2$Downsample*100,"%")
df2$Downsample <- factor(df2$Downsample,levels=unique(df2$Downsample))
df2$Mapped_reads <- NULL

melted <- reshape2::melt(df2)
melted$variable <- ifelse(melted$variable=="Variants","All variants","Somatic variants")

p1=ggplot(melted,aes(x=Downsample,y=value)) +
  theme_classic() +
  geom_point(aes(fill=variable),color='black',shape=21,size=4) +
  # x axis
  theme(axis.text.x=element_text(angle=0, hjust=.5, size = 10)) + 
  theme(axis.title.x = element_text(size=12)) +
  xlab("Downsampling (%)") +
  # y-axis
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  ylab("Number of variants") +
  # scale fill
  scale_fill_manual(values = c("All variants"="grey30","Somatic variants"="grey50")) +
  guides(fill=guide_legend(title="")) +
  # facet
  facet_wrap(~"Mouse #2 - HIP") +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0, face= 'bold')) +
  annotate("text",x="30%",y=27000,label=paste0("p-value = ",pvalue_all),size=3,
           color="grey30",hjust=.7) +
  annotate("text",x="30%",y=24500,label=paste0("p-value = ",pvalue_smt),size=3,
           color="grey50",hjust=.7)
###############################################################################################






###############################################################################################
### Panel B - plot as a coverage

# Calc coverage - downsampling file
downsample$Sample <- "WT_3_HIP"
merged1 <- merge(downsample,mapped,by="Sample")
merged1$coverage <- merged1$Mapped_reads * (merged1$avg_read_len / genome_l)
merged1 <- merged1[,c("Sample","Tissue","Rep","coverage","Variants","Somatic")]


# melt and plot
melted2 <- reshape2::melt(merged1,id.vars = c("Sample","Tissue","Rep","coverage"))
melted2$variable <- ifelse(melted2$variable=="Variants","All variants","Somatic variants")


p2=ggplot(melted2,aes(x=coverage,y=value)) +
  theme_classic() +
  geom_point(aes(fill=variable),color='black',shape=21,size=4) +
  # x axis
  theme(axis.text.x=element_text(angle=0, hjust=.5, size = 10)) + 
  theme(axis.title.x = element_text(size=12)) +
  scale_x_continuous(breaks=seq(0,60, 10)) +
  coord_cartesian(xlim = c(0,65)) +
  xlab("Coverage (X)") +
  # y-axis
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  ylab("Number of variants") +
  # scale fill
  scale_fill_manual(values = c("All variants"="grey30","Somatic variants"="grey50")) +
  guides(fill=guide_legend(title="")) +
  # facet
  facet_wrap(~"Mouse #2 - HIP") +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0, face= 'bold'))
###############################################################################################






###############################################################################################
### Panel C - plot regression line from B and plot values from other samples


# add mapping stat
tt <- merge(tot_variants,mapped,by.x="sample",by.y="Sample")

# calc coverage
tt$coverage <- tt$mapped_reads * (tt$avg_read_len / genome_l)
tt <- tt[,c("sample","Tissue","Rep","coverage","Variants","Somatic")]

# melt
melted_tt <- reshape2::melt(tt,id.vars = c("sample","Tissue","Rep","coverage"))
melted_tt$variable <- ifelse(melted_tt$variable=="Variants","All variants","Somatic variants")
melted_tt$Rep <- paste0("Rep_",melted_tt$Rep-1)

p3=ggplot() + 
  theme_classic() +
  # add regression lines based on previous data
  geom_smooth(data=melted2,aes(x=coverage,y=value,color=variable),method=lm, se=FALSE) +
  # add points from 
  geom_point(data=melted_tt,aes(x=coverage,y=value,fill=Tissue,shape=Rep),size=2,color='black') +
  scale_shape_manual(values = c("Rep_1"=22,"Rep_2"=21,"Rep_3"=23,
                                "Rep_4"=24,"Rep_5"=25)) +
  # x axis
  theme(axis.text.x=element_text(angle=0, hjust=.5, size = 10)) + 
  theme(axis.title.x = element_text(size=12)) +
  scale_x_continuous(breaks=seq(0,60, 10)) +
  coord_cartesian(xlim = c(0,65)) +
  xlab("Coverage (X)") +
  # y-axis
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  ylab("Number of variants") +
  # scale fill
  scale_color_manual(values = c("All variants"="grey30","Somatic variants"="grey50")) +
  guides(color=guide_legend(title="Regression lines")) +
  guides(shape = guide_legend(ncol = 2)) +
  guides(fill = guide_legend(ncol = 2,override.aes = list(shape=21)))
  



riga1=ggarrange(p1,p2,ncol=2,common.legend = T,legend = 'bottom',labels = c("A","B"))
riga2=ggarrange(p3,NULL,ncol=2,labels = "C",widths = c(.9,.1))
pp=ggarrange(riga1,riga2,nrow=2)
ggsave(filename='FigS2.pdf',plot=pp,width=230,height=160,units='mm')

