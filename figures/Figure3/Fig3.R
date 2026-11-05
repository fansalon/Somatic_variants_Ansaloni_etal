#######################################################################################################################################
### Figure 3
#######################################################################################################################################
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)
library(stringr)
library(ggh4x)
library(scales)
library(ggridges)
library(dplyr)
#######################################################################################################################################
source("../../analyses/functions/functions_nanopore.R")
#######################################################################################################################################
### Load input files
# Somatic anno
obj <- readRDS("../../analyses/02_get_somatic/somatic.Rds")
comparison <- obj$tissue_specificity
comparison$Type <- gsub("_vs_","-vs-",comparison$Type)
# Somatic anno+features
bed <- fread("../../analyses/03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
mapped <- fread("../data/mapping_stat.txt",header=T,sep='\t',data.table=F)
# GO results
go_res <- fread("../../analyses/05_GO_tissuespecific_BND_new_version/GO_res_Somatic_tissue_private_BND.txt",header=T,sep='\t',data.table = F)
# genome idx
genome <- fread("../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr",data.table = F)
# total simple rep in mm39 genome
simple_rep <- fread("../../zenodo/input_data/simpleRepeat.tot.bed",data.table = F)
simple_rep <- subset(simple_rep,V1 %in% genome$V1)
simple_rep_occ <- sum(simple_rep$V3-simple_rep$V2) / sum(genome$V2) * 100
#######################################################################################################################################


#######################################################################################################################################

### Select tissue-specific somatic variants
tissue_specific <- bed[grepl("Somatic_tissue_",bed$somatic.anno),]
table(tissue_specific$somatic.anno)
tissue_specific <- subset(tissue_specific,somatic.anno=="Somatic_tissue_private")

# Add comparison type
tissue_specific2 <- merge(tissue_specific,comparison,by.x="ID",by.y="SV_ID",all.x=T,sort=F)

# Expand dataframe based on the "Samples" column (which contains a comma-separated list of the samples where the sv is found)
tissue_specific2_expand <- as.data.frame(tidyr::separate_rows(tissue_specific2, Samples, sep=","))

# Discard samples not matching comparison type
# there might be cases where the samples are KID_6,CTX_2 and the comparisons KID-vs-HIP and CTX-vs-HIP
# in such cases, when the df is expanded, we need to select the matching samples and comparisons
tissue_specific2_expand$tissue_sample <- str_split_fixed(tissue_specific2_expand$Samples,"_",2)[,1]
tissue_specific2_expand$tissue_type <- str_split_fixed(tissue_specific2_expand$Type,"-",3)[,1]
tissue_specific2_expand <- subset(tissue_specific2_expand,tissue_sample==tissue_type)

# Define non-redundant comparison type
tissue_specific2_expand$Type2 <- 
  ifelse(tissue_specific2_expand$Type=="HIP-vs-KID" |tissue_specific2_expand$Type=="KID-vs-HIP","HIP vs\nKID",
         ifelse(tissue_specific2_expand$Type=="CTX-vs-KID" |tissue_specific2_expand$Type=="KID-vs-CTX","CTX vs\nKID",
                ifelse(tissue_specific2_expand$Type=="HIP-vs-CTX" |tissue_specific2_expand$Type=="CTX-vs-HIP","CTX vs\nHIP","error")))

# Modify names in a nicer way ...
tissue_specific2_expand$somatic.anno <- ifelse(tissue_specific2_expand$somatic.anno=="Somatic_tissue_private","Somatic private",
                                               ifelse(tissue_specific2_expand$somatic.anno=="Somatic_tissue_recurrent","Somatic recurrent","error"))
#######################################################################################################################################


#######################################################################################################################################

#### Panel A: norm. number of variants

# Count number of somatic, tissue-specific, variants in each mouse
sv_counts <- as.data.frame(table(tissue_specific2_expand[,c("Samples","Type","somatic.anno","Type2")]))
sv_counts <- subset(sv_counts,Freq>0)

# Normalise by mapped bases
mapped$Samples <- paste0(mapped$Tissue,"_",mapped$Rep)
sv_counts.norm <- merge(sv_counts,mapped,by="Samples",all=T)
sv_counts.norm$Freq.norm <- sv_counts.norm$Freq / sv_counts.norm$mapped_bases * 1^6

# Reorder
df_order <- data.frame("Type2"=c("HIP vs\nKID","CTX vs\nKID","CTX vs\nHIP"))
sv_counts.norm <- merge(df_order,sv_counts.norm,by="Type2",all=T,sort=F)
sv_counts.norm$Type2 <- factor(sv_counts.norm$Type2,levels = unique(sv_counts.norm$Type2))


# T test
counter <- 0
for (somatic in unique(sv_counts.norm$somatic.anno)) {
  d1 <- subset(sv_counts.norm,somatic.anno==somatic)
  for (tp2 in unique(sv_counts.norm$Type2)) {
    counter <- counter + 1
    d2 <- subset(d1,Type2==tp2)
    
    # get counts for the two tissues
    tx1 <- subset(d2,Tissue==unique(d2$Tissue)[1])
    tx2 <- subset(d2,Tissue==unique(d2$Tissue)[2])
    
    # order the same way to make paired t test
    rownames(tx1) <- tx1$Rep
    rownames(tx2) <- tx2$Rep
    tx2 <- tx2[rownames(tx1),]
    
    # t test
    ttest_res <- t.test(tx1$Freq.norm,tx2$Freq.norm,paired=T)
    # save to df
    res <- data.frame("somatic.anno"=somatic,"Type2"=tp2,"y"=max(c(tx1$Freq.norm,tx2$Freq.norm)),"pvalue"=ttest_res$p.value)
    
    # bind
    if (counter==1) {
      rest_stat <- res
    } else {
      rest_stat <- rbind(rest_stat,res)
    }
  }
}
rest_stat$fdr <- p.adjust(rest_stat$pvalue,method = 'fdr')
rest_stat$sign <- ifelse(rest_stat$fdr<0.05&rest_stat$fdr>=0.01,"*",
                         ifelse(rest_stat$fdr<0.01&rest_stat$fdr>=0.001,"**",
                                ifelse(rest_stat$fdr<0.001,"***",NA)))
rest_stat$Type2 <- factor(rest_stat$Type2,levels = unique(rest_stat$Type2))

# Aggregate counts
aggr_counts <- aggregate(sv_counts.norm$Freq,by=list(sv_counts.norm$Type2,sv_counts.norm$Type,sv_counts.norm$somatic.anno),sum)
names(aggr_counts) <- c("Type2","Type","somatic.anno","Tot")
aggr_counts$Freq.norm <- aggregate(sv_counts.norm$Freq.norm,by=list(sv_counts.norm$Type2,sv_counts.norm$Type,sv_counts.norm$somatic.anno),max)[,"x"]

# To have the legend ordered
sv_counts.norm$Type <- factor(sv_counts.norm$Type,levels = unique(sv_counts.norm$Type))

p1=ggplot(subset(sv_counts.norm,somatic.anno=="Somatic private"),aes(x=Type,y=(Freq.norm))) +
  theme_classic() +
  geom_boxplot(aes(fill=Type),alpha=.5,width=.5) +
  # facet
  facet_nested(somatic.anno~Type2,scales = "free") +
  #facet_nested(somatic.anno~Type2,scales = "free") +
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 8, colour = "black", angle = -90, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("CTX-vs-KID"="#619CFF","KID-vs-CTX"="lightblue",
                               "HIP-vs-KID"="#00BA38","KID-vs-HIP"="lightgreen",
                               "CTX-vs-HIP"="#C77CFF","HIP-vs-CTX"="plum1")) +
  scale_x_discrete(labels=c("HIP-vs-CTX"="HIP","KID-vs-CTX"="KID","KID-vs-HIP"="KID",
                            "CTX-vs-HIP"="CTX","CTX-vs-KID"="CTX","HIP-vs-KID"="HIP")) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = 8,angle = 45, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("No. of variants per\n1M mapped bases") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # remove legend
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  # legend spacing
  theme(legend.key.size = unit(0.4, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  geom_point(aes(fill=Type),shape=21,position = position_jitter(width = .2,height = 0)) +
  # add stats
  geom_text(data=subset(rest_stat,somatic.anno=="Somatic private"),aes(label=sign,x=1.5,y=y+1.5e-08),size=6) +
  # add totals
  geom_text(data=subset(aggr_counts,somatic.anno=="Somatic private"),aes(label=Tot,y=Freq.norm+8e-09),size=2) +
  coord_cartesian(ylim = c(2e-08,7.5e-08))
#######################################################################################################################################


#######################################################################################################################################

### Panel B SVTYPE

# Count number of somatic, tissue-specific, variants in each mouse
svt <- as.data.frame(table(tissue_specific2_expand[,c("SVTYPE","Type","somatic.anno","Type2")]))
svt <- subset(svt,Freq>0)

# Calc % for each class
svt$id <- paste0(svt$Type,":",svt$somatic.anno,":",svt$Type2)
svt <- svt[,c("id","SVTYPE","Freq")]
svt <- calcPerc2(svt,group="id",counts="Freq")
svt$Type <- str_split_fixed(svt$id,":",3)[,1]
svt$somatic.anno <- str_split_fixed(svt$id,":",3)[,2]
svt$Type2 <- str_split_fixed(svt$id,":",3)[,3]

# Reorder
svt <- merge(df_order,svt,by="Type2",all=T,sort=F)
svt$Type2 <- factor(svt$Type2,levels = unique(svt$Type2))

# Define who is the first and who is the second bar - for plotting purposes
svt$bar_order <- ifelse(svt$Type=="HIP-vs-KID"|svt$Type=="CTX-vs-KID"|svt$Type=="CTX-vs-HIP","first","second")

# To have the legend ordered
svt$Type <- factor(svt$Type,levels = unique(svt$Type))


p2=ggplot(svt,aes(x=SVTYPE,y=Perc,fill=Type)) +
  theme_classic() +
  geom_bar(stat="identity",position = position_dodge(),color='black',alpha=.7) +
  # facet
  facet_nested(somatic.anno~Type2) +
  #facet_nested(somatic.anno~Type2,scales = "free") +
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 8, colour = "black", angle = -90, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("CTX-vs-KID"="#619CFF","KID-vs-CTX"="lightblue",
                               "HIP-vs-KID"="#00BA38","KID-vs-HIP"="lightgreen",
                               "CTX-vs-HIP"="#C77CFF","HIP-vs-CTX"="plum1")) +
  scale_x_discrete(labels=c("HIP-vs-CTX"="HIP","KID-vs-CTX"="KID","KID-vs-HIP"="KID",
                            "CTX-vs-HIP"="CTX","CTX-vs-KID"="CTX","HIP-vs-KID"="HIP")) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = 8,angle = 45, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% somatic variants") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
#  scale_y_continuous(breaks = seq(0,100,25)) 
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # Add line at 100%
#  geom_hline(yintercept = 100,linetype='dashed',size=.3) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  # Somatic private (BND only)
  geom_text(data=subset(svt, bar_order=="first" & somatic.anno=="Somatic private" & SVTYPE=="BND"),
            aes(label=Freq,y=120),position = position_dodge(1),size=2.5,hjust=.3) +
  geom_text(data=subset(svt, bar_order=="second" & somatic.anno=="Somatic private" & SVTYPE=="BND"),
            aes(label=Freq,y=108),position = position_dodge(1),size=2.5,hjust=-.1) +
  coord_cartesian(ylim = c(0,120)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_hline(yintercept = 100,linetype='dashed',col="grey60")
#######################################################################################################################################


#######################################################################################################################################

### Panel C GO
# Calc log2 enrichment over universe
go_res$log2Enrich <- log2( go_res$Enrichment )

# Add label
go_res$Type <- gsub("_","-",go_res$Type)
go_res$Type2 <- ifelse(go_res$Type=="HIP-vs-KID"|go_res$Type=="KID-vs-HIP","HIP-vs-KID",
                       ifelse(go_res$Type=="CTX-vs-KID"|go_res$Type=="KID-vs-CTX","CTX-vs-KID",
                              ifelse(go_res$Type=="CTX-vs-HIP"|go_res$Type=="HIP-vs-CTX","CTX-vs-HIP","ERR")))

# Reorder
df_orderT <- data.frame("Type2"=c("HIP-vs-KID","CTX-vs-KID","CTX-vs-HIP"))
go_res <- merge(df_orderT,go_res,by="Type2",all=T,sort=F)
go_res$Type2 <- factor(go_res$Type2,levels = unique(go_res$Type2))
df_order2 <- data.frame("Type"=c("HIP-vs-KID","KID-vs-HIP",
                                 "CTX-vs-KID","KID-vs-CTX",
                                 "CTX-vs-HIP","HIP-vs-CTX"))
go_res <- merge(df_order2,go_res,by="Type",all=T,sort=F)
go_res$Type <- factor(go_res$Type,levels = unique(go_res$Type))

# Select the top10 for each category
sel_top <- go_res[,c("term","fdr","features")]
sel_top <- aggregate(sel_top$fdr,by=list(sel_top$term,sel_top$features),mean)
names(sel_top) <- c("term","features","fdr")
top10 <- as.data.frame(sel_top %>%
                         group_by(features) %>%
                         slice_min(fdr, n = 10,with_ties=F))
top10$id <- paste0(top10$term,":",top10$features)


# Get top10
go_res$id <- paste0(go_res$term,":",go_res$features)
dt <- subset(go_res,id %in% top10$id)

# plot only terms which are fdr<0.05 in at least one comparison
to_keep <- unique(subset(dt,fdr<0.05)[,"term"])
dt2 <- subset(dt,term %in% to_keep)

# Plot
dt2$term <- ifelse(dt2$term=="nucleobase-containing compound metabolic process",
                       "nucleob.-contain. compound metab. proc.",dt2$term)
dt2$term <- ifelse(dt2$term=="negative regulation of nucleobase-containing compound metabolic process",
                   "neg. reg. of nucleob.-contain. compound metab. proc.",dt2$term)
#dt2$term <- gsub("negative regulation","neg. reg.",dt2$term)
p3=ggplot(dt2,aes(x=term,y=-log10(fdr))) +
  theme_bw() +
  # white background of the facet
  theme(strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2))) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.20, 0.15))) +
  # to increase readibility
  # geom_vline(xintercept = 1:nrow(dt),linetype='dashed',linewidth=.1) +
  # line at fdr 0.05
  geom_point(aes(size=log2Enrich,fill=log10(numDEInCat)),shape=21) +
  geom_hline(yintercept = -log10(0.05),color='red',linetype='dashed') +
  scale_fill_gradient(low = "blue",high = "red") +
  facet_nested(features~Type2+Type,scales = "free_y",space = "free_y") +
  theme(strip.text.x = element_text(size = 7.3, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 7.5, colour = "black", angle = -90, face= 'bold')) +
  # y axis
  xlab("none") +
  theme(axis.text.y=element_text(size = 7,angle = 0, hjust = 1)) +
  theme(axis.title.y = element_text(size=0)) +
  # x axis
  ylab("FDR (-log10)") +
  theme(axis.text.x=element_text(size = 8,angle = 0, hjust = 1)) +
  theme(axis.title.x = element_text(size=10)) +
  # legend title
  guides(size=guide_legend(title="Enrichment\n(log2)")) +
  guides(fill=guide_colourbar(title="No. of genes\n(log10)")) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.001, 'cm')) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=9)) + 
  scale_size(range = c(1.5, 6), name = "Enrichment (log2)")
#######################################################################################################################################


#######################################################################################################################################

### Panel D: overlap simple rep

# Select BND
tissue_specific2_expand_bnd <- subset(tissue_specific2_expand,SVTYPE=="BND")

# count breakpoint in SR, and not SR if only 1 breakpoint is in SR
tissue_specific2_expand_bnd$SR_count <- ifelse(is.na(tissue_specific2_expand_bnd$SR),
                                               "No SR/No SR",tissue_specific2_expand_bnd$SR)
sr_class_1 <- data.frame(SR_count = str_split_fixed(tissue_specific2_expand_bnd$SR_count,"/",2)[,1],
                         Type     = tissue_specific2_expand_bnd$Type,
                         Type2    = tissue_specific2_expand_bnd$Type2)
sr_class_2 <- data.frame(SR_count = str_split_fixed(tissue_specific2_expand_bnd$SR_count,"/",2)[,2],
                         Type     = tissue_specific2_expand_bnd$Type,
                         Type2    = tissue_specific2_expand_bnd$Type2)
sr_class <- rbind(sr_class_1,sr_class_2)

# convert to binary
sr_class$SR_count <- ifelse(sr_class$SR_count!="No SR" & grepl(")n",sr_class$SR_count),"SR",sr_class$SR_count)

# count
sr_ct <- as.data.frame(table(sr_class[,c("SR_count","Type","Type2")]))
sr_ct <- subset(sr_ct,Freq>0)
sr_ct <- calcPerc2(sr_ct,"Type","Freq")
sr_ct <- subset(sr_ct,SR_count=="SR")

# Get "control" from mm39
mm39 <- data.frame("Type"="mm\n39",
                   "SR_count"="SR",
                   "Type2"="mm\n39",
                   "Freq"=sum(simple_rep$V3-simple_rep$V2),
                   "Tot"=sum(genome$V2) ,
                   "Perc"=simple_rep_occ)

# bind
tot <- rbind(sr_ct,mm39)

# Reorder
tot <- merge(df_order,tot,by="Type2",all=T,sort=F)
tot$Type2 <-factor(tot$Type2,levels=unique(tot$Type2))

# calculate stats
mm39 <- subset(tot,Type=="mm\n39")
all <- subset(tot,Type!="mm\n39")
counter <- 0
for (tp in unique(all$Type)) {
  counter <- counter + 1
  d2 <- subset(all,Type==tp)
  res <- prop.test( c(d2$Freq,unique(mm39$Freq)), c(d2$Tot,unique(mm39$Tot)), alternative = "greater")
  stat_res <- data.frame("Type"=tp,"Type2"=unique(d2$Type2),"pvalue"=res$p.value)
  # bind
  if (counter==1) {
    stat_res2 <- stat_res
  } else {
    stat_res2 <- rbind(stat_res2,stat_res)
  }
}
stat_res2$fdr <- p.adjust(stat_res2$pvalue,method = 'fdr')
stat_res2$sign <- ifelse(stat_res2$fdr<0.05&stat_res2$fdr>=0.01,"*",
                         ifelse(stat_res2$fdr<0.01&stat_res2$fdr>=0.001,"**",
                                ifelse(stat_res2$fdr<0.001,"***",NA)))

p4=ggplot(tot,aes(x=Type,y=Perc,fill=Type)) +
  theme_classic() +
  geom_bar(stat="identity",color='black',width=.7) +
  # facet
  #facet_nested(~somatic.anno+Type2,scales = "free") +
  facet_nested("Somatic private\nBND"~Type2,scales = "free",space = "free") +
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 8, colour = "black", angle = -90, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("CTX-vs-KID"="#619CFF","KID-vs-CTX"="lightblue",
                               "HIP-vs-KID"="#00BA38","KID-vs-HIP"="lightgreen",
                               "CTX-vs-HIP"="#C77CFF","HIP-vs-CTX"="plum1")) +
  scale_x_discrete(labels=c("HIP-vs-CTX"="HIP","KID-vs-CTX"="KID","KID-vs-HIP"="KID",
                            "CTX-vs-HIP"="CTX","CTX-vs-KID"="CTX","HIP-vs-KID"="HIP")) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = 8,angle = 45, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% in simple repeats") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  #scale_y_continuous(breaks = seq(0,100,25)) +
  coord_cartesian(ylim = c(0,5)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # Add line at 100%
  #  geom_hline(yintercept = 100,linetype='dashed',size=.3) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  guides(fill = FALSE) +
  # scale alpha
  scale_alpha_manual(values = c("High complexity"=.2,"Low complexity"=.8)) +
  # add stat
  geom_text(data=stat_res2,aes(label=sign,y=3),size=4,alpha=1) 
#  # add label
#  geom_text(data=subset(tot,Type!="mm\n39"),aes(label=Freq,y=Perc+1),alpha=1,size=2.4)
#######################################################################################################################################


#######################################################################################################################################



riga1=ggarrange(p1,p2,p4,ncol=3,labels = c("A","B","D"),common.legend = T,legend = 'bottom')
riga2=ggarrange(p3,labels = "C")

ff=ggarrange(riga1,riga2,nrow=2,ncol=1,heights = c(.33,.66))
ggsave(filename='Fig3.pdf',plot=ff,width=210,height=200,units='mm')


