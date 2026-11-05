#######################################################################################################################################
### Figure 4
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
# Somatic anno+features
bed <- fread("../../analyses/03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
mapped <- fread("../data/mapping_stat.txt",header=T,sep='\t',data.table=F)
# Simple rep occupancy in mm39
simple_rep_mm39 <- fread("../../zenodo/input_data/simpleRepeat.tot.bed",data.table = F)
genome <- fread("../../zenodo/input_data/GRCm39.vM28.primary_assembly.genome.fa.fai.noscaff.chr",data.table = F)
# SR shuffling res (motif len)
sr_overlap <- fread("../../analyses/09_common_DELs_SR/DEL_somatic_common_SR_shuffling_1000.txt",header=T,sep='\t',data.table = F)
# SR shuffling res (motif type)
rep_c <- fread("../../analyses/10_common_DELs_SR_N/DEL_somatic_common_SRmotif_shuffling_1000.txt",data.table=F,sep='\t')
rep_c2 <- fread("../../analyses/11_common_DELs_SR_2N/DEL_somatic_common_SRmotif_shuffling_1000.txt",data.table=F,sep='\t')
# GO results
go_res <- fread("../../analyses/08_GO_common_DELs_new_version/GO_res_Somatic_common_DEL.txt",header=T,sep='\t',data.table = F)
#######################################################################################################################################


#######################################################################################################################################

### Select common somatic SVs
common <- bed[grepl("Somatic_common_",bed$somatic.anno),]

# Count in how many mice the variant is found
common$no_mice <- str_count(common$Samples,",") + 1

# Modify names in a nicer way ...
common$somatic.anno <- ifelse(common$somatic.anno=="Somatic_common_private","Somatic private",
                              ifelse(common$somatic.anno=="Somatic_common_recurrent","Somatic recurrent",
                                     ifelse(common$somatic.anno=="Somatic_common_unclassified","Somatic unclassified","error")))
#######################################################################################################################################



#######################################################################################################################################

#### Panel A: number of SVs

# Count number of somatic, tissue-specific, SVs in each mouse
sv_counts <- as.data.frame(table(common[,c("no_mice","somatic.anno")]))
sv_counts <- subset(sv_counts,Freq>0)

# order
sv_counts <- sv_counts[order(sv_counts$no_mice),]

# Modify label
sv_counts$no_mice <- paste0(sv_counts$no_mice,"/5")
sv_counts$no_mice <- factor(sv_counts$no_mice,levels = unique(sv_counts$no_mice))

p1=ggplot(sv_counts,aes(x=no_mice,y=(Freq),fill=somatic.anno)) +
  theme_classic() +
  geom_bar(stat="identity",color='black',width=.8,alpha=.7) +
  # scale colors
  scale_fill_manual(values = c("Somatic private"="lightpink","Somatic unclassified"="mistyrose",
                               "Somatic recurrent"="firebrick")) +
  # x axis
  xlab("No. of mice") +
  theme(axis.text.x=element_text(size = 8,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=10)) +
  # y axis
  ylab("No. of loci hosting\nsomatic common variant") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # labels
  geom_text(aes(label=Freq,y=Freq+400),size=3) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

#######################################################################################################################################


#######################################################################################################################################

### Panel B SVTYPE

common2 <- subset(common,somatic.anno=="Somatic private"|somatic.anno=="Somatic recurrent")

# Count number of somatic, tissue-specific, SVs in each mouse
svt <- as.data.frame(table(common2[,c("SVTYPE","somatic.anno")]))
svt <- subset(svt,Freq>0)

# Calc % for each class
svt <- calcPerc2(svt,group="somatic.anno",counts="Freq")

# modify name
svt$somatic.anno <- gsub(" ","\n",svt$somatic.anno)


p2=ggplot(svt,aes(x=SVTYPE,y=Perc,fill=somatic.anno)) +
  theme_classic() +
  geom_bar(stat="identity",position = position_dodge(),color='black',alpha=.7) +
  # facet
  facet_nested(~somatic.anno) +
  #facet_nested(somatic.anno~Type2,scales = "free") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 9, colour = "black", angle = -90, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("Somatic\nprivate"="lightpink","Somatic\nrecurrent"="firebrick")) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = 7.5,angle = 30, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("% of loci hosting\nsomatic common variant") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # Add line at 100%
  geom_hline(yintercept = 100,linetype='dashed',size=.3) +
  # legend position
  theme(legend.position = "none") +
  # Somatic private (BND only)
  geom_text(data=subset(svt, SVTYPE=="DEL"|SVTYPE=="INS"),
            aes(label=Freq,y=Perc+7),position = position_dodge(1),size=3,hjust=.6)
#######################################################################################################################################


#######################################################################################################################################

#### Panel C: SVLEN
del <- subset(common2,SVTYPE=="DEL")
ins <- subset(common2,SVTYPE=="INS")
tt=rbind(del,ins)

# annotate
mm <- aggregate(tt$SVLEN,by=list(tt$SVTYPE,tt$somatic.anno),median)
names(mm) <- c("SVTYPE","somatic.anno","median")

p3=ggplot(tt, aes(x = SVLEN, y = somatic.anno,fill=somatic.anno)) + 
  geom_density_ridges(scale = .5,bandwidth=.08,alpha=.5,
                      jittered_points = F,point_size=.1,rel_min_height=0.00000000000001) + 
  theme_classic() +
  # x axis
  xlab("SV length [bp]") +
  theme(axis.text.x=element_text(size = 8,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=10)) +
  # facet
  facet_nested(~SVTYPE) +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  # y axis
  ylab("no title") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=0)) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  theme(legend.title = element_blank()) +
  # scale colors
  scale_fill_manual(values = c("Somatic private"="lightpink","Somatic recurrent"="firebrick")) +
  # legend
  theme(legend.position = "none") +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  #  geom_vline(xintercept = 25,linetype='dashed',color="black") +
  # annotate("text",x=80,y=0.8,label="25 bp",size=3) +
  geom_text(data=subset(mm,somatic.anno=="Somatic recurrent"),aes(label=paste0("median = ",median," bp"),x=800,y=2.4),size=2.8) +
  geom_text(data=subset(mm,somatic.anno!="Somatic recurrent"),aes(label=paste0("median = ",median," bp"),x=800,y=1.3),size=2.8)
#######################################################################################################################################

#######################################################################################################################################

### Panel D - GO
# Calc log2 enrichment over universe
go_res$log2Enrich <- log2( go_res$Enrichment )

# Select the top10 for each category
sel_top <- go_res[,c("term","fdr","features")]
sel_top <- aggregate(sel_top$fdr,by=list(sel_top$term,sel_top$features),mean)
names(sel_top) <- c("term","features","fdr")
top10 <- as.data.frame(sel_top %>%
                         group_by(features) %>% 
                         slice_min(fdr, n = 15,with_ties=F))
top10$id <- paste0(top10$term,":",top10$features)


# Get top10
go_res$id <- paste0(go_res$term,":",go_res$features)
dt <- subset(go_res,id %in% top10$id)


# plot only terms which are fdr<0.05 in at least one comparison
to_keep <- unique(subset(dt,fdr<0.1)[,"term"])
dt2 <- subset(dt,term %in% to_keep)

# Prep for plot
dt2$Type <- ifelse(dt2$somatic.anno=="Somatic_common_private","Somatic\nprivate",
                  ifelse(dt2$somatic.anno=="Somatic_common_recurrent","Somatic\nrecurrent","ERR"))

# order
dt2      <- dt2[order(dt2$fdr,decreasing = F),]
dt2$term <- factor(dt2$term,levels = unique(dt2$term))

# Plot
p4=ggplot(dt2,aes(x=term,y=-log10(fdr))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2))) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 1, 2),
                     expand = expansion(mult = c(0.20, 0.15))) +
  # line at fdr 0.05
  geom_point(aes(size=log2Enrich,fill=log10(numDEInCat)),shape=21) +
  geom_hline(yintercept = -log10(0.1),color='red',linetype='dashed') +
  scale_fill_gradient(low = "blue",high = "red") +
  facet_nested(features~Type,scales = "free_y",space = "free_y") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 9, colour = "black", angle = -90, face= 'bold')) +
  # y axis
  xlab("none") +
  theme(axis.text.y=element_text(size = 9,angle = 0, hjust = 1)) +
  theme(axis.title.y = element_text(size=0)) +
  # x axis
  ylab("FDR (-log10)") +
  theme(axis.text.x=element_text(size = 8,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=10)) +
  # legend title
  guides(size=guide_legend(title="Enrichment\n(log2)")) +
  guides(fill=guide_colourbar(title="No. of genes\n(log10)")) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.001, 'cm')) + 
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=9))
#######################################################################################################################################

#######################################################################################################################################

### Panel E: overlap simple rep
tot <- as.data.frame(table(tt[,c("overlaps.simplerep","somatic.anno","SVTYPE")]))
tot$id <- paste0(tot$somatic.anno,"@",tot$SVTYPE)
tot <- calcPerc2(tot,"id","Freq")
tot$id <- NULL


### Get "control" from mm39
sr_occu <- sum(simple_rep_mm39$V3-simple_rep_mm39$V2)
genome_len <- sum(genome$V2)
control <- data.frame(overlaps.simplerep= c("TRUE","FALSE"),
                      somatic.anno      = c("mm39","mm39"),
                      SVTYPE            = c("mm39","mm39"),
                      Freq              = c(sr_occu,genome_len-sr_occu),
                      Tot               = genome_len,
                      Perc              = c(sr_occu/genome_len*100,(genome_len-sr_occu)/genome_len*100))

# bind
tot <- rbind(tot,control)

# rename
tot$somatic.anno <- ifelse(tot$somatic.anno=="Somatic private","Private",
                           ifelse(tot$somatic.anno=="Somatic recurrent","Recurrent",as.character(tot$somatic.anno)))
tot$somatic.anno <- factor(tot$somatic.anno,levels=unique(tot$somatic.anno))


p5=ggplot(subset(tot,SVTYPE!="INS"),aes(x=somatic.anno,y=Perc,fill=somatic.anno,alpha=overlaps.simplerep)) +
  theme_classic() +
  geom_bar(stat="identity",color='black',width=.5) +
  # facet
  facet_nested(~"DEL",scales = "free",space = "free") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("Private"="lightpink","Recurrent"="firebrick")) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = 8,angle = 30, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab(" % of loci hosting\nsomatic common DEL") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  scale_y_continuous(breaks = seq(0,100,25)) +
  coord_cartesian(ylim = c(0,110)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  theme(legend.title = element_text(size=10)) +
  guides(fill = FALSE) +
  guides(alpha=guide_legend(title="Simple\nrepeats",nrow=2,byrow=TRUE)) +
  # scale alpha
  scale_alpha_manual(values = c("FALSE"=.1,"TRUE"=.8)) +
  # add stat
  geom_text(data=subset(tot,overlaps.simplerep=="TRUE" &SVTYPE=="DEL"),aes(label="***",y=Perc-15),alpha=1,size=4) +
  # add label
  geom_text(data=subset(tot,overlaps.simplerep=="TRUE" &SVTYPE=="DEL"),aes(label=Freq,y=Perc-5),alpha=1,size=2.8) +
  geom_text(data=subset(tot,overlaps.simplerep=="FALSE" &SVTYPE=="DEL"),aes(label=Freq,y=100+6),alpha=1,size=2.8)
#######################################################################################################################################

#######################################################################################################################################


### Panel F: SR motif length
# get % expected
sr_overlap$Exp_perc <- sr_overlap$Exp / sr_overlap$Tot * 100
# get % sd
sr_overlap$sd_perc <- sr_overlap$sd / sr_overlap$Tot * 100

# classify >=10 as others
sr_overlap$SR_len_class <- ifelse(sr_overlap$SR_motif_len==1,"1",
                                  ifelse(sr_overlap$SR_motif_len==2,"2",
                                         ifelse(sr_overlap$SR_motif_len==3,"3",
                                                ifelse(sr_overlap$SR_motif_len==4,"4",
                                                       ifelse(sr_overlap$SR_motif_len==5,"5",
                                                              ifelse(sr_overlap$SR_motif_len==6,"6",
                                                                     ifelse(sr_overlap$SR_motif_len==7,"7",
                                                                            ifelse(sr_overlap$SR_motif_len==8,"8",
                                                                                   ifelse(sr_overlap$SR_motif_len==9,"9","Others")))))))))
# Aggregate counts according to SR_len_class
sr_overlap2 <- aggreCounts2(sr_overlap)

# Define significant
sr_overlap2$sign <- ifelse(sr_overlap2$fdr<0.05&sr_overlap2$fdr>=0.01,"*",
                           ifelse(sr_overlap2$fdr<0.01&sr_overlap2$fdr>=0.001,"**",
                                  ifelse(sr_overlap2$fdr<0.001,"***",NA)))
# Consider as significant only if Obs % > 5%
sr_overlap2$sign <- ifelse(sr_overlap2$Obs_perc>5,sr_overlap2$sign,NA)


# Prepare to plot
obs <- sr_overlap2[,c("somatic.anno","SR_len_class","Obs_perc","Obs","fdr","Enrichment","sign")]
names(obs)[names(obs)=="Obs_perc"] <- "Perc"
names(obs)[names(obs)=="Obs"] <- "Counts"
obs$sd <- NA
obs$shuf <- "Observed"
obs$sign <- ifelse(obs$Enrichment>0,obs$sign,NA)
exp <- sr_overlap2[,c("somatic.anno","SR_len_class","Exp_perc","Exp","fdr","Enrichment","sign","sd")]
names(exp)[names(exp)=="Exp_perc"] <- "Perc"
names(exp)[names(exp)=="Exp"] <- "Counts"
exp$sign <- ifelse(exp$Enrichment<0,exp$sign,NA)
exp$shuf <- "Expected"
tot <- rbind(obs,exp)
tot$shuf <- factor(tot$shuf,levels = unique(tot$shuf))

# Reorder 1
df_order2 <- data.frame("SR_len_class"=c("1","2","3","4","5","6","7","8","9","Others"))
tot <- merge(df_order2,tot,by="SR_len_class",all=T,sort=F)
tot$SR_len_class <- factor(tot$SR_len_class,levels = unique(tot$SR_len_class))

# add totals
tt <- aggregate(tot$Counts,by=list(tot$somatic.anno,tot$shuf),sum)
names(tt) <- c("somatic.anno","shuf","n")

p6=ggplot(tot,aes(x=SR_len_class,y=Counts,fill=somatic.anno,alpha=shuf)) + 
  geom_bar(stat="identity",color='black',position = position_dodge()) +
  geom_errorbar(aes(ymin=Counts-sd, ymax=Counts+sd),
                width=0.35,position = position_dodge(1),linewidth=.8) +
  scale_alpha_manual(values = c("Observed"=1,"Expected"=.4)) +
  theme_classic() +
  facet_nested(~somatic.anno,scales = "free",space = "free") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 9, colour = "black", angle = -90, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("Somatic private"="lightpink","Somatic recurrent"="firebrick")) +
  # x axis
  xlab("Simple repeat motif length") +
  theme(axis.text.x=element_text(size = 8,angle = 90, hjust = 1,vjust=.5)) +
  theme(axis.title.x = element_text(size=10)) +
  # y axis
  ylab("No. loci hosting common somatic\nDELs in simple repeats") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # legend pos and title
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  # Add line at 100%
  #  geom_hline(yintercept = 100,linetype='dashed',size=.3) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  guides(fill = FALSE) +
  geom_text(aes(label=sign,y=Counts+50),size=4,alpha=1) +
  # add totals
  geom_text(data=subset(tt,shuf=="Observed"),aes(x=8,y=,5500,label=paste0("Tot = ",n)),size=3,alpha=1)
#######################################################################################################################################

#######################################################################################################################################


### Panel G: STR classification in (N)n and (NN)n


########### (N)n
# rename some cols
names(rep_c)[names(rep_c)=="Perc"] <- "Obs_perc"
# Exp perc
rep_c$Exp_perc <- rep_c$Exp / rep_c$Tot * 100
# sd perc
rep_c$sd_perc <- rep_c$sd / rep_c$Tot * 100
# define stat signif.
rep_c$sign <- ifelse(rep_c$fdr<0.05&rep_c$fdr>=0.01,"*",
                     ifelse(rep_c$fdr<0.01&rep_c$fdr>=0.001,"**",
                            ifelse(rep_c$fdr<0.001,"***",NA)))
# Consider as significant only if Obs % > 5%
rep_c$sign <- ifelse(rep_c$Obs_perc>5,rep_c$sign,NA)

# Prepare to plot
obs <- rep_c[,c("somatic.anno","SR","Obs_perc","Obs","fdr","Enrichment","sign")]
names(obs)[names(obs)=="Obs_perc"] <- "Perc"
names(obs)[names(obs)=="Obs"] <- "Counts"
obs$sd <- NA
obs$shuf <- "Observed"
obs$sign <- ifelse(obs$Enrichment<0,NA,obs$sign)
exp <- rep_c[,c("somatic.anno","SR","Exp_perc","Exp","fdr","Enrichment","sign","sd")]
names(exp)[names(exp)=="Exp_perc"] <- "Perc"
names(exp)[names(exp)=="Exp"] <- "Counts"
exp$shuf <- "Expected"
exp$sign <- ifelse(exp$Enrichment>0,NA,exp$sign)
tot <- rbind(obs,exp)
tot$shuf <- factor(tot$shuf,levels = unique(tot$shuf))

# add label
tot$Type2 <- "(N)n"

########### (NN)n
# rename some cols
names(rep_c2)[names(rep_c2)=="Perc"] <- "Obs_perc"
# Exp perc
rep_c2$Exp_perc <- rep_c2$Exp / rep_c2$Tot * 100
# sd perc
rep_c2$sd_perc <- rep_c2$sd / rep_c2$Tot * 100
# define stat signif.
rep_c2$sign <- ifelse(rep_c2$fdr<0.05&rep_c2$fdr>=0.01,"*",
                      ifelse(rep_c2$fdr<0.01&rep_c2$fdr>=0.001,"**",
                             ifelse(rep_c2$fdr<0.001,"***",NA)))
# Consider as significant only if Obs % > 5%
rep_c2$sign <- ifelse(rep_c2$Obs_perc>5,rep_c2$sign,NA)


# Prepare to plot
obs <- rep_c2[,c("somatic.anno","SR","Obs_perc","Obs","fdr","Enrichment","sign")]
names(obs)[names(obs)=="Obs_perc"] <- "Perc"
names(obs)[names(obs)=="Obs"] <- "Counts"
obs$sd <- NA
obs$shuf <- "Observed"
obs$sign <- ifelse(obs$Enrichment<0,NA,obs$sign)
exp <- rep_c2[,c("somatic.anno","SR","Exp_perc","Exp","fdr","Enrichment","sign","sd")]
names(exp)[names(exp)=="Exp_perc"] <- "Perc"
names(exp)[names(exp)=="Exp"] <- "Counts"
exp$shuf <- "Expected"
exp$sign <- ifelse(exp$Enrichment>0,NA,exp$sign)
tot2 <- rbind(obs,exp)
tot2$shuf <- factor(tot2$shuf,levels = unique(tot2$shuf))

# add label
tot2$Type2 <- "(NN)n"

totf <- rbind(tot,tot2)

p7=ggplot(subset(totf,Perc>1),aes(x=SR,y=Counts,fill=somatic.anno,alpha=shuf)) + 
  geom_bar(stat="identity",color='black',position = position_dodge()) +
  geom_errorbar(aes(ymin=Counts-sd, ymax=Counts+sd),
                width=0.35,position = position_dodge(1),linewidth=.8) +
  scale_alpha_manual(values = c("Observed"=1,"Expected"=.4)) +
  theme_classic() +
  facet_nested(Type2~somatic.anno,scales = "free_x",independent = "x") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 9, colour = "black", angle = -90, face= 'bold')) +
  # scale colors
  scale_fill_manual(values = c("Somatic private"="lightpink","Somatic recurrent"="firebrick")) +
  # x axis
  xlab("Simple repeat motif") +
  theme(axis.text.x=element_text(size = 8,angle = 90, hjust = 1,vjust=.5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("No. loci hosting common somatic\nDELs in simple repeats") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=10)) +
  # legend pos and title
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  # Add line at 100%
  #  geom_hline(yintercept = 100,linetype='dashed',size=.3) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  guides(fill = FALSE) +
  geom_text(aes(label=sign,y=Counts+20),size=4,alpha=1) +
  coord_cartesian(ylim = c(0,1800))
#######################################################################################################################################

#######################################################################################################################################


### Panel H - length of SR containing DELs
common_SR_DEL <- subset(common2,final.anno=="simple_rep" & SVTYPE =="DEL")

# Get SR length
common_SR_DEL$SR_end <- as.numeric(str_split_fixed(common_SR_DEL$SR_ID,"-",3)[,3])
common_SR_DEL$SR_start <- as.numeric(str_split_fixed(common_SR_DEL$SR_ID,"-",3)[,2])
common_SR_DEL$SR_len <- common_SR_DEL$SR_end - common_SR_DEL$SR_start 

# median
med_lab <- aggregate(common_SR_DEL$SR_len,by=list(common_SR_DEL$somatic.anno),median)
names(med_lab) <- c("somatic.anno","median")

# change names
common_SR_DEL$somatic.anno <- ifelse(common_SR_DEL$somatic.anno=="Somatic private","Somatic\nprivate",
                                     ifelse(common_SR_DEL$somatic.anno=="Somatic recurrent","Somatic\nrecurrent","ERR"))

p8=ggplot(common_SR_DEL, aes(x = SR_len, y = somatic.anno,fill=somatic.anno)) + 
  geom_density_ridges(scale = .5,bandwidth=.08,alpha=.5,
                      jittered_points = F,point_size=.1,rel_min_height=0.00000000000001) + 
  theme_classic() +
  # x axis
  xlab("SR length [bp]") +
  theme(axis.text.x=element_text(size = 8,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=10)) +
  # facet
  facet_nested(~"Simple repeats\ncontaining DELs") +
  theme(strip.text.x = element_text(size = 9, colour = "black", angle = 0, face= 'bold')) +
  # y axis
  ylab("no title") +
  theme(axis.text.y=element_text(size = 8)) +
  theme(axis.title.y = element_text(size=0)) +
  # legend spacing
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.spacing.x = unit(.1, 'cm'),legend.spacing.y = unit(.1, 'cm')) +
  theme(legend.title = element_blank()) +
  # scale colors
  scale_fill_manual(values = c("Somatic\nprivate"="lightpink","Somatic\nrecurrent"="firebrick")) +
  # legend
  theme(legend.position = "none") +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  geom_text(data=subset(med_lab,somatic.anno=="Somatic recurrent"),aes(label=paste0("median=",median," bp"),x=2500,y=2.4),size=2.8) +
  geom_text(data=subset(med_lab,somatic.anno!="Somatic recurrent"),aes(label=paste0("median=",median," bp"),x=2500,y=1.3),size=2.8)


#######################################################################################################################################

#######################################################################################################################################








#######################################################################################################################################

#######################################################################################################################################
p11=ggarrange(p1,nrow=2,heights = c(.985,.015))
p22=ggarrange(p2,nrow=2,heights = c(.975,.025))
p77=ggarrange(p7,nrow=2,heights = c(.943,.057))
p67=ggarrange(p6,p77,ncol=2,labels = c("F","G"),common.legend = T,legend = 'bottom',widths = c(.48,.52))
p88=ggarrange(p8,nrow=2,heights = c(.855,.145),labels = "H")

riga1=ggarrange(p11,p22,p3,ncol=3,widths = c(.8,.8,1.4),labels = c("A","B","C"))
riga2=ggarrange(p4,p5,ncol=2,widths = c(.7,.3),labels = c("D","E"))
riga3=ggarrange(p67,p88,ncol=2,widths = c(.78,.22))
ff=ggarrange(riga1,riga2,riga3,nrow=3,ncol=1,heights = c(.8,1.1,1.1))
ggsave(filename='Fig4.pdf',plot=ff,width=210,height=222,units='mm')
#######################################################################################################################################



