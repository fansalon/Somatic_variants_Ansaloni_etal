#######################################################################################################################################
### Figure S5
#######################################################################################################################################
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(data.table)
library(scales)
library(stringr)
library(reshape2)
library(ggh4x)
#######################################################################################################################################
source("../../analyses/functions/functions_nanopore.R")
#######################################################################################################################################
### Graphical parameters
textbarsize=2.5
axistexty=8
axistitley=12
#######################################################################################################################################
mei <- fread("../../analyses/13_tldr/input/tldr_res.txt",data.table=F,sep='\t',header=T)
mei_f <- fread("../../analyses/13_tldr/tldr_res_filtered.txt.bed",data.table=F,sep='\t',header=T)
mei_c <- fread("../../analyses/13_tldr/tldr_res_filtered_reads_summary.txt",data.table=F,sep='\t',header=T)
genomicf <- fread("../../analyses/14_annotate_MEI/res_observed_expected.txt",data.table=F,sep='\t',header=T)
mei_gf <- fread("../../analyses/14_annotate_MEI/genomic_feature_anno.txt",data.table=F,sep='\t',header=T)
mapped <- fread("../data/mapping_stat.txt",data.table=F,sep='\t',header=T)
#######################################################################################################################################


###### Panel A

### Count different type of MEIs

count_mei <- data.frame(c(
             "Tot"                      = nrow(mei),
             "PASS"                     = nrow(subset(mei,Filter=="PASS")),
             "PASS_TSD"                 = nrow(subset(mei,Filter=="PASS" & !is.na(mei$TSD))),
             "PASS_TSD_nonref"          = nrow(subset(mei,Filter=="PASS" & !is.na(mei$TSD) & is.na(mei$NonRef) )),
             "PASS_TSD_nonref_somatic"  = nrow(subset(mei_f,is.somatic))
))
colnames(count_mei) <- "Count"
count_mei$x <- rownames(count_mei)


# plot
count_mei$x <- factor(count_mei$x,levels = unique(count_mei$x))
p1=ggplot(count_mei,aes(x=x,y=Count)) +
  theme_classic() +
  geom_bar(stat='identity',color='black',aes(fill=x),alpha=.8) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = axistexty,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # modify x axis labels
  scale_x_discrete(labels=c("Tot"                     = "Tot",
                            "PASS"                    = "PASS",
                            "PASS_TSD"                = "PASS\nTSD",
                            "PASS_TSD_nonref"         = "PASS\nTSD\nunknown",
                            "PASS_TSD_nonref_somatic" = "PASS\nTSD\nunknown\nsomatic")) +
  # y axis
  ylab("No. of loci hosting MEI") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=axistitley)) +
  #scale_y_continuous(trans='log2', labels = trans_format("log2", math_format(2^.x))) +
  scale_y_continuous(trans='log2') +
  # label
  geom_text(aes(label=Count,y=Count+(Count*50/100)),size=3) +
  # legend
  guides(fill="none") +
  # scale colors
  scale_fill_manual( values = c("Tot"                      = "grey20",
                                "PASS"                     = "grey50",
                                "PASS_TSD"                 = "grey80",
                                "PASS_TSD_nonref"          = "grey100",
                                "PASS_TSD_nonref_somatic"  = hue_pal()(3)[3]))
#######################################################################################################################################

#######################################################################################################################################



###### Panel B

somatic_mei <- subset(mei_f,is.somatic==TRUE)

### Type of TE
count_te <- as.data.frame(table(somatic_mei$TE))
count_te$perc <- count_te$Freq / sum(count_te$Freq) * 100
count_te$Var1 <- str_split_fixed(count_te$Var1,":",2)[,2]

# label
count_te$lbl <- ifelse(count_te$Freq>=7,count_te$Freq,NA)

# modify IAPLTR1a_MM to IAPLTR1a
count_te$Var1 <- ifelse(count_te$Var1=="IAPLTR1a_MM","IAPLTR1a",count_te$Var1)

p2=ggplot(count_te,aes(x="Somatic MEI",y=Freq,fill=Var1)) +
  theme_classic() +
  geom_bar(stat='identity',color='black',alpha=.8,width=.6) +
  # label
  geom_text(aes(label=lbl,x=1,y=Freq), position = position_stack(vjust = .5),size=3) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = axistexty,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("No. of loci hosting\nsomatic MEI") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=axistitley)) +
  # scale colors
  scale_fill_manual( values = c("L1MdA_I"       = "steelblue1",
                                "L1MdGf_I"      = "dodgerblue2",
                                "L1MdTf_I"      = "royalblue2",
                                "L1MdTf_II"     = "slateblue2",
                                "L1MdTf_III"    = "purple1",
                                "IAPLTR1a"      = "green3",
                                "B1"            = "firebrick1",
                                "B2"            = "indianred")
  )  +
  # legend
  guides(fill=guide_legend(title = '',nrow=4,byrow=TRUE)) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=8)) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm'))
#######################################################################################################################################

#######################################################################################################################################



###### Panel C

# genomic features
genomicf$sd_perc <- genomicf$sd / genomicf$Tot

# aggregate <1
genomicf$annotation <- ifelse(genomicf$Obs==0,"Other",genomicf$annotation)

# separate obs and exp
obs <- genomicf[,c("annotation","Obs","Obs_perc","zscore","fdr")]
obs$sd_perc <- NA
names(obs)[names(obs)=="Obs"] <- "count"
names(obs)[names(obs)=="Obs_perc"] <- "perc"
obs$Type <- "Somatic MEI"

exp <- genomicf[,c("annotation","Exp","Exp_perc","zscore","fdr","sd_perc")]
names(exp)[names(exp)=="Exp"] <- "count"
names(exp)[names(exp)=="Exp_perc"] <- "perc"
exp$Type <- "Exp"

# get tot
tot <- rbind(obs,exp)
tot$Type <- factor(tot$Type,levels = unique(tot$Type))

# get label
tot$lbl <- paste0( round(tot$count), " (",round(tot$perc),"%)")
tot$lbl <- ifelse(tot$annotation=="Other",NA,tot$lbl)

# sign
tot$sign <- ifelse(tot$fdr<0.05 & tot$fdr>=0.01,"*",
                  ifelse(tot$fdr<0.01 & tot$fdr>=0.001,"**",
                         ifelse(tot$fdr<0.001,"***",NA)))
tot$sign <- ifelse(tot$Type!="Exp",tot$sign,NA)

p3=ggplot(tot,aes(x=Type,y=perc,fill=annotation)) +
  theme_classic() +
  geom_bar(stat='identity',color='black',alpha=.6) +
  # label
  geom_text(aes(label=lbl), position = position_stack(vjust = .5),size=3) +
  # label sign
  geom_text(aes(label=sign), position = position_stack(vjust = .6),size=5) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = axistexty,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% of loci hosting\nsomatic MEI") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=axistitley)) +
  # legend
  guides(fill=guide_legend(title = '',nrow=4,byrow=TRUE)) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=8)) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  # scale colors
  scale_fill_manual( values = c("Distal Intergenic" = "lightsalmon",
                                "Intron"            = "lightgoldenrod1",
                                "Promoter"          = "lightyellow1",
                                "Other"             = "black"))
#######################################################################################################################################

#######################################################################################################################################


###### Panel D
# ALL L1 somatic MEI

# add info on ALT to somatic MEI
somatic_mei2 <- merge(somatic_mei,mei_c[,c("ID",colnames(mei_c)[grepl("_ALT",colnames(mei_c))])],by="ID")

##########################################################################################
#### 1
#### Prepare variables

# Separate common/tissue-spec from recurrent/private etc
somatic_mei2$Type2 <- str_split_fixed(somatic_mei2$Type,"_",2)[,1]
somatic_mei2$Type <- str_split_fixed(somatic_mei2$Type,"_",2)[,2]

# order
somatic_mei2 <- somatic_mei2[order(somatic_mei2$Type),]
somatic_mei2$ID <- factor(somatic_mei2$ID,levels = unique(somatic_mei2$ID))

# MEI length to kb
somatic_mei2$LengthIns <- somatic_mei2$LengthIns/1000

# MAF in %
somatic_mei2$mean_MAF <- somatic_mei2$mean_MAF*100


##########################################################################################
#### 2
#### Get new values

# calc mean number of supporting reads
somatic_mei2$mean_ALT <- apply(somatic_mei2[,grepl("_ALT",colnames(somatic_mei2))],1,sum) / somatic_mei2$NumSamples

# calc number of tissues MEIs are found in each mouse
for (j in 2:6) {
  somatic_mei2[,paste0("Mouse ",j-1)] <- rowSums(somatic_mei2[,grepl("_ALT",colnames(somatic_mei2)) & 
                                                                grepl(paste0("WT_",j),colnames(somatic_mei2))]>=1)
}

# calc number of mice where MEIs are found in each tissue
tissues <- c("CTX","KID","HIP")
for (tx in tissues) {
  somatic_mei2[,tx] <- rowSums(somatic_mei2[,grepl("_ALT",colnames(somatic_mei2)) & grepl(tx,colnames(somatic_mei2))]>=1)
}

# add genomic features
somatic_mei2 <- merge(somatic_mei2,mei_gf[,c("ID","annotation","gene_name","distanceTx","distanceTSS")],by='ID',all.x=T)

# There are MEI that are not tissuespec or common. They are found in all the three tissues
# but never in the same mice
somatic_mei2$Type2 <- ifelse(somatic_mei2$Type2=="Tissue-specific" & somatic_mei2$CTX>0 & somatic_mei2$KID>0 & somatic_mei2$HIP>0,
                             "Unclassified",somatic_mei2$Type2)
somatic_mei2$Type <- ifelse(somatic_mei2$Type2=="Unclassified","unclassified",somatic_mei2$Type)
##########################################################################################
#### 3
#### Select columns of interest
somatic_mei22 <- somatic_mei2[,c("ID","TE","LengthIns","mean_ALT","mean_MAF","NumSamples",
                                 "Mouse 1","Mouse 2","Mouse 3","Mouse 4","Mouse 5",
                                 "CTX","HIP","KID","annotation","Type","Type2","full_length")]

# convert annotation to numeric
somatic_mei22$annotation2 <- ifelse(somatic_mei22$annotation=="Promoter",3,
                                    ifelse(somatic_mei22$annotation=="Intron",3,
                                           ifelse(somatic_mei22$annotation=="Distal Intergenic",3,"ERR")))
somatic_mei22$annotation2 <- as.numeric(somatic_mei22$annotation2)


##########################################################################################
#### 4
#### melt

# melt
somatic_mei222 <- melt(somatic_mei22)


##########################################################################################
#### 5
#### fill

somatic_mei222$fill <- "no"

# annotation2, for anno; full_length for LengthIns
somatic_mei222$fill <- ifelse(somatic_mei222$variable=="annotation2",somatic_mei222$annotation,somatic_mei222$fill)
somatic_mei222$full_length <- ifelse(somatic_mei222$full_length==TRUE,"full-length","non-full-length")
somatic_mei222$fill <- ifelse(somatic_mei222$variable=="LengthIns",somatic_mei222$full_length,somatic_mei222$fill)

# type, for NumSamples, Mouse and tissues
somatic_mei222$fill <- ifelse(somatic_mei222$variable=="NumSamples"|somatic_mei222$variable=="Mouse 1"|somatic_mei222$variable=="Mouse 2"|
                                somatic_mei222$variable=="Mouse 3"|somatic_mei222$variable=="Mouse 4"|somatic_mei222$variable=="Mouse 5"|
                                somatic_mei222$variable=="CTX"|somatic_mei222$variable=="HIP"|somatic_mei222$variable=="KID",somatic_mei222$Type,somatic_mei222$fill)


##########################################################################################
### 6 plot
###

# Only L1
somatic_mei222_l1 <- somatic_mei222[grepl("L1",somatic_mei222$TE),]

# Rename facet
somatic_mei222_l1$variable <- ifelse(somatic_mei222_l1$variable=="LengthIns","MEI len\n(kb)",
                                     ifelse(somatic_mei222_l1$variable=="mean_ALT","Supp.\nreads",
                                            ifelse(somatic_mei222_l1$variable=="mean_MAF","MAF\n(%)",
                                                   ifelse(somatic_mei222_l1$variable=="NumSamples","Numb.\nsamples",
                                                          ifelse(somatic_mei222_l1$variable=="annotation2","Anno-\ntation",
                                                                 ifelse(somatic_mei222_l1$variable=="Mouse 1","Mouse1\n# tissues",
                                                                        ifelse(somatic_mei222_l1$variable=="Mouse 2","Mouse2\n# tissues",
                                                                               ifelse(somatic_mei222_l1$variable=="Mouse 3","Mouse3\n# tissues",
                                                                                      ifelse(somatic_mei222_l1$variable=="Mouse 4","Mouse4\n# tissues",
                                                                                             ifelse(somatic_mei222_l1$variable=="Mouse 5","Mouse5\n# tissues",
                                                                                                    ifelse(somatic_mei222_l1$variable=="KID","KID\n# mice",
                                                                                                           ifelse(somatic_mei222_l1$variable=="CTX","CTX\n# mice",
                                                                                                                  ifelse(somatic_mei222_l1$variable=="HIP","HIP\n# mice",
                                                                                                                         as.character(somatic_mei222_l1$variable))))))))))))))

# factor
somatic_mei222_l1$variable <- factor(somatic_mei222_l1$variable,levels = unique(somatic_mei222_l1$variable))

# add gene name
gene_anno <- somatic_mei2[somatic_mei2$ID %in% somatic_mei222_l1$ID,c("ID","annotation","gene_name","distanceTx")]
gene_anno$Gene <- ifelse(gene_anno$distanceTx<10000,gene_anno$gene_name,as.character(gene_anno$ID))
somatic_mei222_l1 <- merge(somatic_mei222_l1,gene_anno[,c("ID","Gene")],by="ID",all.x=T,sort=F)

somatic_mei222_l1$fill <- ifelse(somatic_mei222_l1$fill=="unclassified","unclassif.",somatic_mei222_l1$fill)

p4=ggplot(somatic_mei222_l1,aes(x=Gene,y=value,fill=fill)) +
  theme_bw() +
  geom_bar(stat="identity",color='black') +
  facet_nested(Type2~variable,scales = "free",space = "free_y") +
  theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 9, colour = "black", angle = -90, face= 'bold')) +
  coord_flip() +
  # x axis
  xlab("somatic LINE-1 insertions") +
  theme(axis.text.x=element_text(size = axistexty,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("none") +
  theme(axis.text.y=element_text(size = 5)) +
  theme(axis.title.y = element_text(size=axistitley)) +
  # legend
  guides(fill=guide_legend(title = '',nrow=1,byrow=TRUE)) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  # scale fill
  scale_fill_manual( values = c("full-length"="indianred1",
                                "non-full-length"="seagreen3",
                                "no"="lightblue",
                                "recurrent"="black","unclassif."="grey50","private"="white",
                                "Intron"="orange","Promoter"="purple1","Distal Intergenic"="lightyellow"))
#######################################################################################################################################

#######################################################################################################################################

p11=ggarrange(p1,nrow=2,heights = c(.86,.14))
riga1=ggarrange(p11,p2,p3,nrow=1,widths = c(.46,.27,.27),labels = c("A","B","C"))
riga2=ggarrange(p4,nrow=1,labels = c("D"),legend = 'bottom')
ff=ggarrange(riga1,riga2,NULL,nrow=3,ncol=1,heights = c(.3,.7,0))
ggsave(filename='FigS5.pdf',plot=ff,width=210,height=297,units='mm')




