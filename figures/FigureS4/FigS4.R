#######################################################################################################################################
### Figure S4 - distribution genomic features - common
#######################################################################################################################################
rm(list = ls())
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggh4x)
library(dplyr)
#######################################################################################################################################
# genomic feature results
genfeat <- fread("../../analyses/04_annotate_genomic_features/final_anno_table_ALL_1000.txt",header=T,sep='\t',data.table = F)
genfeat_SR <- fread("../../analyses/04_annotate_genomic_features/final_anno_table_SimpleRep_1000.txt",header=T,sep='\t',data.table = F)
genfeat_normal <- fread("../../analyses/04_annotate_genomic_features/final_anno_table_normal_1000.txt",header=T,sep='\t',data.table = F)
# Somatic anno
obj <- readRDS("../../analyses/02_get_somatic/somatic.Rds")
# Somatic anno+features
bed <- fread("../../analyses/03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
# GO results
go_res_SR <- fread("../../analyses/08_GO_common_DELs/GO_res_Simple_rep_Somatic_common_DEL.txt",header=T,sep='\t',data.table = F)
go_res_normal <- fread("../../analyses/08_GO_common_DELs/GO_res_normal_Somatic_common_DEL.txt",header=T,sep='\t',data.table = F)
#######################################################################################################################################

# Select genomic features for common variants
genfeat <- subset(genfeat,somatic.anno!="Somatic_tissue_private")
genfeat$name <- "All"
genfeat_SR <- subset(genfeat_SR,somatic.anno!="Somatic_tissue_private")
genfeat_SR$name <- "Simple\nrepeats"
genfeat_normal <- subset(genfeat_normal,somatic.anno!="Somatic_tissue_private")
genfeat_normal$name <- "Non simple\nrepeats"

# Bind
final <- rbind(genfeat,genfeat_SR,genfeat_normal)
final$sd_perc <- final$sd / final$Tot *100

# Get only DELs
final <- subset(final,SVTYPE=="DEL")

# Format for plot
obs <- final[,c("annotation","zscore","Obs_perc","fdr","somatic.anno","name")]
names(obs)[3] <- c("Perc")
obs$sd_perc <- NA
obs$Type <- "Observed"
exp <- final[,c("annotation","zscore","Exp_perc","fdr","somatic.anno","name","sd_perc")]
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

# Get pval ***
combined$sign <- ifelse(combined$fdr<0.05&combined$fdr>=0.01,"*",
                        ifelse(combined$fdr<0.01&combined$fdr>=0.001,"**",
                               ifelse(combined$fdr<0.001,"***",NA)))
#####################################################################################################################
# Sum SVs for plotting totals
genfeat$name <- "All"
genfeat_SR$name <- "Simple\nrepeats"
genfeat_normal$name <- "Non simple\nrepeats"
tt <- rbind(genfeat,genfeat_normal,genfeat_SR)
tt <- subset(tt,SVTYPE=="DEL")
totals <- aggregate(tt$Obs,by=list(tt$somatic.anno,tt$name),sum)
names(totals) <- c("somatic.anno","name","Tot")
#####################################################################################################################

# rename
combined$somatic.anno <- ifelse(combined$somatic.anno=="Somatic_common_private","Somatic private (DEL)",
                                ifelse(combined$somatic.anno=="Somatic_common_recurrent","Somatic recurrent (DEL)","ERR"))
totals$somatic.anno <- ifelse(totals$somatic.anno=="Somatic_common_private","Somatic private (DEL)",
                                ifelse(totals$somatic.anno=="Somatic_common_recurrent","Somatic recurrent (DEL)","ERR"))

# not significant in Perc<1%
combined$sign <- ifelse(combined$Perc<2,NA,combined$sign)

p1=ggplot(data=subset(combined,name=="All"),aes(x=annotation,y=Perc,fill=somatic.anno,alpha=Type)) +
  theme_classic() +
  geom_bar(stat='identity',position = position_dodge(),color='black') +
  # error bars
  geom_errorbar(aes(ymin=Perc-sd_perc, ymax=Perc+sd_perc),
                width=0.35,position = position_dodge(1),linewidth=.8) +
  # facet
  facet_nested(name~somatic.anno) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 7, colour = "black", face= 'bold')) +
  scale_alpha_manual(values = c("Observed"=.8,"Expected"=.3)) +
  # scale colors
  scale_fill_manual(values = c("Somatic private (DEL)"="lightpink","Somatic recurrent (DEL)"="firebrick")) +
  # x axis
  xlab("none") +
  #theme(axis.text.x=element_text(size = 8,angle = 90, hjust = 1,vjust = .5)) +
  theme(axis.text.x=element_text(size = 8,angle = 45, hjust = 1,vjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% loci hosting somatic common DEL") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  guides(fill = "none") + 
  # add stat
  geom_text(aes(label=sign,y=Perc+2),size=5,alpha=1) +
  # add label
  geom_text(data=subset(totals,name=="All"),aes(x=1.5,y=60,label=paste0("Tot=",Tot)),alpha=1,size=3) +
  coord_cartesian(ylim = c(0,60))
#######################################################################################################################################
#######################################################################################################################################


### Panel B - divided in overlappingg or not overlapping simpe rep
p2=ggplot(data=subset(combined,name!="All"),aes(x=annotation,y=Perc,fill=somatic.anno,alpha=Type)) +
  theme_classic() +
  geom_bar(stat='identity',position = position_dodge(),color='black') +
  # error bars
  geom_errorbar(aes(ymin=Perc-sd_perc, ymax=Perc+sd_perc),
                width=0.35,position = position_dodge(1),linewidth=.8) +
  # facet
  facet_nested(name~somatic.anno) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 7, colour = "black", face= 'bold')) +
  scale_alpha_manual(values = c("Observed"=.8,"Expected"=.3)) +
  # scale colors
  scale_fill_manual(values = c("Somatic private (DEL)"="lightpink","Somatic recurrent (DEL)"="firebrick")) +
  # x axis
  xlab("none") +
  #theme(axis.text.x=element_text(size = 8,angle = 90, hjust = 1,vjust = .5)) +
  theme(axis.text.x=element_text(size = 8,angle = 45, hjust = 1,vjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% somatic variants (common DEL)") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  guides(fill = "none") + 
  # add stat
  geom_text(aes(label=sign,y=Perc+3),size=5,alpha=1) +
  # add label
  geom_text(data=subset(totals,name!="All"),aes(x=1.5,y=64,label=paste0("Tot=",Tot)),alpha=1,size=3) +
  coord_cartesian(ylim = c(0,65))

#######################################################################################################################################

#######################################################################################################################################


### Panel B: GO simple rep
# Calc log2 enrichment over universe
go_res_SR$log2Enrich <- log2( (go_res_SR$intersection_size/go_res_SR$query_size) / (go_res_SR$term_size/go_res_SR$effective_domain_size) )


# Select the top10 for each category
sel_top <- go_res_SR[,c("term_name","p.adj","annotation","Type")]
sel_top <- aggregate(sel_top$p.adj,by=list(sel_top$term_name,sel_top$annotation,sel_top$Type),mean)
names(sel_top) <- c("term_name","annotation","Type","p.adj")
top10_rec <- as.data.frame(subset(sel_top,Type=="Somatic_common_recurrent") %>%
                             group_by(annotation) %>% 
                             slice_min(p.adj, n = 15,with_ties=F))
top10_rec$id <- paste0(top10_rec$term_name,":",top10_rec$annotation)
top10_priv <- as.data.frame(subset(sel_top,Type!="Somatic_common_recurrent") %>%
                             group_by(annotation) %>% 
                             slice_min(p.adj, n = 15,with_ties=F))
top10_priv$id <- paste0(top10_priv$term_name,":",top10_priv$annotation)

top10 <- unique(rbind(top10_rec,top10_priv))


# Plot
go_res_SR$id <- paste0(go_res_SR$term_name,":",go_res_SR$annotation)
dt <- subset(go_res_SR,id %in% top10$id)
dt$Type <- ifelse(dt$Type=="Somatic_common_private","Somatic private",
                  ifelse(dt$Type=="Somatic_common_recurrent","Somatic recurrent","ERR"))
pBOH=ggplot(dt,aes(x=term_name,y=-log10(p.adj))) +
  theme_classic() +
  coord_flip() +
  # to increase readibility
  geom_vline(xintercept = 1:nrow(dt),linetype='dashed',linewidth=.1) +
  # line at fdr 0.05
  geom_hline(yintercept = -log10(0.05),color='red',linetype='dashed') +
  geom_point(aes(size=log2Enrich,fill=log10(intersection_size)),shape=21) +
  scale_fill_gradient(low = "blue",high = "red") +
  facet_nested(annotation~Type,scales = "free_y",space = "free_y") +
  theme(strip.text.x = element_text(size = 7.5, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 7.5, colour = "black", angle = -90, face= 'bold')) +
  # y axis
  xlab("none") +
  theme(axis.text.y=element_text(size = 9,angle = 0, hjust = 1)) +
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
  ggtitle("Simple repeats") + 
  theme(plot.title = element_text(size=12,hjust = 0.5,face = 'bold'))
#######################################################################################################################################


#######################################################################################################################################

ggsave(filename='FigS4.pdf',plot=p1,width=210,height=198/2,units='mm')

