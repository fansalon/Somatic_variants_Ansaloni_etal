#######################################################################################################################################
#######################################################################################################################################
rm(list = ls())
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggh4x)
library(dplyr)
#######################################################################################################################################
# TE enrichment results
te_enrich <- fread("../../analyses/12_tissuespec_BND_TEs/BND_overlapTEs_res_boot1000.txt",header=T,sep='\t',data.table = F)
# genomic feature results
genfeat <- fread("../../analyses/04_annotate_genomic_features/final_anno_table_ALL_1000.txt",header=T,sep='\t',data.table = F)
genfeat_SR <- fread("../../analyses/04_annotate_genomic_features/final_anno_table_SimpleRep_1000.txt",header=T,sep='\t',data.table = F)
genfeat_normal <- fread("../../analyses/04_annotate_genomic_features/final_anno_table_normal_1000.txt",header=T,sep='\t',data.table = F)
#######################################################################################################################################

### Panel A: TE enrichment
te_enrich$Type2 <- ifelse(te_enrich$comparison=="CTX-vs-HIP"|te_enrich$comparison=="HIP-vs-CTX","CTX-vs-HIP",
                          ifelse(te_enrich$comparison=="CTX-vs-KID"|te_enrich$comparison=="KID-vs-CTX","CTX-vs-KID",
                                 ifelse(te_enrich$comparison=="HIP-vs-KID"|te_enrich$comparison=="KID-vs-HIP","HIP-vs-KID","ERR")))

# select top
topXX <- unique(subset(te_enrich,Obs>15)[,"TE"])
selected <- te_enrich[te_enrich$TE %in% topXX,]


# fill
selected$fill <- ifelse(selected$Enrichment>0,"Enriched","Depleted")

# Reorder cols
df_order <- data.frame("Type2"=c("HIP-vs-KID","CTX-vs-KID","CTX-vs-HIP"))
selected <- merge(df_order,selected,by="Type2",all=T,sort=F)
selected$Type2 <- factor(selected$Type2,levels = unique(selected$Type2))
df_order2 <- data.frame("Type"=c("HIP-vs-KID","KID-vs-HIP",
                                 "CTX-vs-KID","KID-vs-CTX",
                                 "CTX-vs-HIP","HIP-vs-CTX"))
selected <- merge(df_order2,selected,by.x="Type",by.y="comparison",all=T,sort=F)
selected$Type <- factor(selected$Type,levels = unique(selected$Type))

# fill based on fdr
selected$Significant <- ifelse(selected$fdr<0.05 & selected$Obs_perc>5,"TRUE","FALSE")

# order by Obs_perc
selected <- selected[order(selected$Obs_perc,decreasing = F),]
selected$TE <- factor(selected$TE,levels = unique(selected$TE))

#p0=ggplot(selected,aes(x=TE,y=-log10(fdr))) +
p0=ggplot(selected,aes(x=TE,y=zscore)) +
  theme_classic() +
  coord_flip(ylim = c(-4,4)) +
  # to increase readibility
  geom_vline(xintercept = 1:length(unique(selected$TE)),linetype='dashed',linewidth=.1) +
  # line at fdr 0.05
  geom_hline(yintercept = 0,color='red',linetype='dashed') +
  geom_point(aes(size=Obs_perc,fill=Significant),shape=21,alpha=.8) +
  scale_fill_manual(values = c("TRUE"="red","FALSE"="grey")) +
 # scale_fill_gradient(low = "blue",high = "red") +
  facet_nested(~Type2+Type,scales = "free_y",space = "free_y") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 9, colour = "black", angle = -90, face= 'bold')) +
  # y axis
  xlab("BND breakpoints annotated over TEs") +
  theme(axis.text.y=element_text(size = 8,angle = 0, hjust = 1)) +
  theme(axis.title.y = element_text(size=12)) +
  # x axis
  ylab("zscore") +
  theme(axis.text.x=element_text(size = 8,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=12)) +
  # legend title
  guides(size=guide_legend(title="% BNDs")) +
  #guides(fill=guide_legend(title="")) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.001, 'cm')) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=9))





#######################################################################################################################################

### Panel B: genomif  feat
# Select genomic features for tissue specific BND
genfeat <- subset(genfeat,somatic.anno=="Somatic_tissue_private")
genfeat_SR <- subset(genfeat_SR,somatic.anno=="Somatic_tissue_private")
genfeat_normal <- subset(genfeat_normal,somatic.anno=="Somatic_tissue_private")


### Aggregate annotations found in <2%
inp_ls <- list(genfeat,genfeat_SR,genfeat_normal)
names(inp_ls) <- c("All","Simple\nrepeats","Non simple\nrepeats")
out_ls <- list()
i <- 0
for (inp in inp_ls) {
  i <- i + 1
  # <2%  are classified as Other
  inp$annotation <- ifelse(inp$Obs_perc<=2,"Other",inp$annotation)
  # Obs
  genfeat2 <- aggregate(inp$Obs_perc,by=list(inp$annotation,inp$type),sum)
  names(genfeat2) <- c("annotation","type","Obs_perc")
  # Exp
  genfeat2$Exp_perc <- aggregate(inp$Exp_perc,by=list(inp$annotation,inp$type),sum)[,"x"]
  # Sd
  # Get sd in %
  inp$sd_perc <- inp$sd / inp$Tot * 100
  # Get var (var is cum, sd is not)
  inp$var_perc <- inp$sd_perc^2
  genfeat2$var_perc <- aggregate(inp$var_perc,by=list(inp$annotation,inp$type),sum)[,"x"]
  # convert to sd
  genfeat2$sd_perc <- sqrt(genfeat2$var_perc)
  genfeat2$var_perc <- NULL
  # zscore
  genfeat2$zscore <- (genfeat2$Obs_perc - genfeat2$Exp_perc) / genfeat2$sd_perc
  # pvalue
  genfeat2$pvalue <- 2*pnorm(-abs(genfeat2$zscore))
  # fdr 
  genfeat2$fdr <- p.adjust(genfeat2$pvalue,method = 'fdr')
  
  # bind
  genfeat2$name <- names(inp_ls)[[i]]
  if (i==1){
    final <- genfeat2
  } else {
    final <- rbind(final,genfeat2)
  }
}



### Format for plot
obs <- final[,c("annotation","zscore","Obs_perc","fdr","type","name")]
names(obs)[3] <- c("Perc")
obs$sd_perc <- NA
obs$Type <- "Observed"
exp <- final[,c("annotation","zscore","Exp_perc","fdr","type","name","sd_perc")]
names(exp)[3] <- c("Perc")
exp$Type <- "Expected"
combined <- rbind(obs,exp)

# Modify name
combined$annotation <- ifelse(combined$annotation=="Distal Intergenic","Intergenic",combined$annotation)

# factor
combined$Type <- factor(combined$Type,levels = unique(combined$Type))

# assign fdr to the highest bar
combined$fdr <- ifelse(combined$zscore>0 & combined$Type=="Expected",NA,combined$fdr)
combined$fdr <- ifelse(combined$zscore<0 & combined$Type=="Observed",NA,combined$fdr)

# Get facet
combined$type <- gsub("_","-",combined$type)
combined$Type2 <- ifelse(combined$type=="HIP-vs-KID"|combined$type=="KID-vs-HIP","HIP-vs-KID",
                         ifelse(combined$type=="CTX-vs-KID"|combined$type=="KID-vs-CTX","CTX-vs-KID",
                                ifelse(combined$type=="CTX-vs-HIP"|combined$type=="HIP-vs-CTX","CTX-vs-HIP","ERR")))

# Get pval ***
combined$sign <- ifelse(combined$fdr<0.05&combined$fdr>=0.01,"*",
                        ifelse(combined$fdr<0.01&combined$fdr>=0.001,"**",
                               ifelse(combined$fdr<0.001,"***",NA)))

# Reorder
df_order <- data.frame("Type2"=c("HIP-vs-KID","CTX-vs-KID","CTX-vs-HIP"))
combined <- merge(df_order,combined,by="Type2",all=T,sort=F)
combined$Type2 <- factor(combined$Type2,levels = unique(combined$Type2))
df_order2 <- data.frame("type"=c("HIP-vs-KID","KID-vs-HIP",
                                 "CTX-vs-KID","KID-vs-CTX",
                                 "CTX-vs-HIP","HIP-vs-CTX"))
combined <- merge(df_order2,combined,by="type",all=T,sort=F)
combined$type <- factor(combined$type,levels = unique(combined$type))


# Modify intergenic label (Distal Intergenic --> Intergenic)
combined$annotation <- gsub("Distal ","",combined$annotation)


# Sum SVs for plotting totals
genfeat$name <- "All"
genfeat_SR$name <- "Simple\nrepeats"
genfeat_normal$name <- "Non simple\nrepeats"
tt <- rbind(genfeat,genfeat_normal,genfeat_SR)
totals <- aggregate(tt$Obs,by=list(tt$type,tt$name),sum)
names(totals) <- c("type","name","Tot")
totals$type <- gsub("_","-",totals$type)
totals$Type2 <- ifelse(totals$type=="HIP-vs-KID"|totals$type=="KID-vs-HIP","HIP-vs-KID",
                       ifelse(totals$type=="CTX-vs-KID"|totals$type=="KID-vs-CTX","CTX-vs-KID",
                              ifelse(totals$type=="CTX-vs-HIP"|totals$type=="HIP-vs-CTX","CTX-vs-HIP","ERR")))

# Reorder
totals <- merge(df_order,totals,by="Type2",all=T,sort=F)
totals$Type2 <- factor(totals$Type2,levels = unique(totals$Type2))
totals <- merge(df_order2,totals,by="type",all=T,sort=F)
totals$type <- factor(totals$type,levels = unique(totals$type))


p1=ggplot(data=combined,aes(x=annotation,y=Perc,fill=type,alpha=Type)) +
  theme_classic() +
  geom_bar(stat='identity',position = position_dodge(),color='black') +
  # error bars
  geom_errorbar(aes(ymin=Perc-sd_perc, ymax=Perc+sd_perc),
                width=0.35,position = position_dodge(1),linewidth=.8) +
  # facet
  facet_nested(name~Type2+type) +
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 8, colour = "black", face= 'bold')) +
  scale_alpha_manual(values = c("Observed"=.8,"Expected"=.3)) +
  # x axis
  xlab("none") +
  #theme(axis.text.x=element_text(size = 8,angle = 90, hjust = 1,vjust = .5)) +
  theme(axis.text.x=element_text(size = 10,angle = 90, hjust = 1,vjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% somatic variants\n(tissue-specific private BND)") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # Add line at 100%
  #  geom_hline(yintercept = 100,linetype='dashed',size=.3) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  # add stat
  geom_text(aes(label=sign),size=5) +
  # scale colors
  scale_fill_manual(values = c("CTX-vs-KID"="#619CFF","KID-vs-CTX"="lightblue",
                               "HIP-vs-KID"="#00BA38","KID-vs-HIP"="lightgreen",
                               "CTX-vs-HIP"="#C77CFF","HIP-vs-CTX"="plum1")) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  # add label
  geom_text(data=totals,aes(x=5,y=30,label=paste0("Tot=",Tot)),alpha=1,size=2.5)
#######################################################################################################################################

#######################################################################################################################################

#ff=ggarrange(p0,p1,labels = c("A","B"),nrow=2,heights = c(.61,.39))
ff=ggarrange(p0,p1,labels = c("A","B"),nrow=2,heights = c(.61,.39))

#ggsave(filename='FigS3.pdf',plot=ff,width=210,height=297*2/3,units='mm')
ggsave(filename='FigS3.pdf',plot=ff,width=210,height=297,units='mm')

