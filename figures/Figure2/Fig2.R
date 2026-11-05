#######################################################################################################################################
### Figure 2
#######################################################################################################################################
rm(list = ls())
library(ComplexUpset)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(stringr)
#######################################################################################################################################
### Graphical parameters
# Upset
textbarsize=1.9
axistexty=10
axistitley=12
titlesize=12
set.size_name="Variants per sample"
set.size_xtext=7
set.size_xtitle=9
#######################################################################################################################################
obj <- readRDS("../../analyses/02_get_somatic/somatic.Rds")
bed <- fread("../../analyses/03_annotate_vcf_somatic/parsed_vcf_anno_somatic.bed.gz",header=T,sep='\t',data.table = F)
#######################################################################################################################################


#######################################################################################################################################

#### Panel E - upset 3/3 variants
common2 <- obj$common_somatic

# Prams for labelling the vertical bars
ylabb="Common\nsomatic variants 0/0\n(3/3 tissues)"
base_ann.ls <- list(
  ylabb=(intersection_size(mapping=aes(fill=intersection,color=intersection),
                           text=list(size=textbarsize),bar_number_threshold=1,color='black') ) +
    scale_fill_manual(breaks = c("1-2-3-4-5","2","1-2"),
                      labels=c("Somatic recurrent","Somatic private","Somatic unclassified"),
                      values = c(
                        # somatic private in white
                        "1"="white","2"="white","3"="white","4"="white","5"="white",
                        # somatic recurrent in black
                        "1-2-3-4-5"="black",
                        # others in default grey - this is a fake only to get this label on the legend
                        "1-2"="grey50"
                      )) +
    # guides(fill = "none") +
    theme(legend.title = element_blank()) +
    theme(legend.key.size = unit(3, 'mm')) +
    theme(legend.position = 'top') +
    theme(axis.text.y=element_text(size = axistexty-1)) +
    theme(axis.title.y=element_text(size = axistitley))
  #    ggtitle("CICCI") +
  #    theme(plot.title = element_text(size=titlesize,hjust=0.5,face='bold'))
)
names(base_ann.ls) <- ylabb


panel_A=ComplexUpset::upset(common2,
                            intersect = colnames(common2),
                            name="",
                            labeller=ggplot2::as_labeller(c(
                              "n_mouse_WT2"="Mouse 1 - 3/3 tissues",
                              "n_mouse_WT3"="Mouse 2 - 3/3 tissues",
                              "n_mouse_WT4"="Mouse 3 - 3/3 tissues",
                              "n_mouse_WT5"="Mouse 4 - 3/3 tissues",
                              "n_mouse_WT6"="Mouse 5 - 3/3 tissues"
                          #    "germinal"="0/1 or 1/1 in any tissue or mouse"
                            )),
                            base_annotations=base_ann.ls,
                            set_sizes=FALSE,
                            # set_sizes=upset_set_size(position='right') +
                            #   #scale_y_continuous(trans=reverse_log_trans(),labels=log10) +
                            #   #ylab('log10(set size)')
                            #   ylab(set.size_name) +
                            #   theme(axis.text.x=element_text(size = set.size_xtext)) +
                            #   theme(axis.title.x=element_text(size = set.size_xtitle)),
                            # width_ratio=0.13,
                            min_size=0)
#######################################################################################################################################




#######################################################################################################################################
#### Panel B - upset 1/3 and 2/3 variants
tx_spec2 <- obj$tissue_specific_somatic

# Prams for labelling the vertical bars
ylabb="Tissue-specific\nsomatic variants 0/0\n(1/3 or 2/3 tissues)"
base_ann.ls <- list(
  ylabb=(intersection_size(mapping=aes(fill=intersection,color=intersection),
                           text=list(size=textbarsize),bar_number_threshold=1,color='black') ) +
    scale_fill_manual(breaks = c("1-2-3-4-5","2","1-2"),
                      labels=c("Somatic recurrent","Somatic private","Somatic unclassified"),
                      values = c(
                        # somatic private in white
                        "1"="white","2"="white","3"="white","4"="white","5"="white",
                        # somatic recurrent in black
                        "1-2-3-4-5"="black",
                        # others in default grey - this is a fake only to get this label on the legend
                        "1-2"="grey50"
                      )) +
    guides(fill = "none") +
    theme(legend.title = element_blank()) +
    theme(legend.key.size = unit(3, 'mm')) +
   # theme(legend.position = 'top') +
    theme(axis.text.y=element_text(size = axistexty-1)) +
    theme(axis.title.y=element_text(size = axistitley))
  #    ggtitle("CICCI") +
  #    theme(plot.title = element_text(size=titlesize,hjust=0.5,face='bold'))
)
names(base_ann.ls) <- ylabb


panel_B=ComplexUpset::upset(tx_spec2,
                            intersect = colnames(tx_spec2),
                            name="",
                            labeller=ggplot2::as_labeller(c(
                              "Mouse_2"="Mouse 1 - 1/3 or 2/3 tissues",
                              "Mouse_3"="Mouse 2 - 1/3 or 2/3 tissues",
                              "Mouse_4"="Mouse 3 - 1/3 or 2/3 tissues",
                              "Mouse_5"="Mouse 4 - 1/3 or 2/3 tissues",
                              "Mouse_6"="Mouse 5 - 1/3 or 2/3 tissues"
                              #   "germinal"="0/1 or 1/1 in any tissue or mouse"
                            )),
                            base_annotations=base_ann.ls,
                            # set_sizes=upset_set_size(position='right') +
                            #   #scale_y_continuous(trans=reverse_log_trans(),labels=log10) +
                            #   #ylab('log10(set size)')
                            #   ylab(set.size_name) +
                            #   theme(axis.text.x=element_text(size = set.size_xtext-3)) +
                            #   theme(axis.title.x=element_text(size = set.size_xtitle)),
                            # width_ratio=0.13,
                            # guides='over',
                            set_sizes=FALSE,
                            min_size=0
)





#######################################################################################################################################




#######################################################################################################################################

#### Panel C supporting reads

# Get somatic only
somatic <- bed[!is.na(bed$somatic.anno) & !grepl("unclassified",bed$somatic.anno),]
table(somatic$somatic.anno)

# Get supporting reads only
somatic_supp <- somatic[,grepl("somatic.anno",colnames(somatic)) |
                          grepl("_ALT",colnames(somatic)) &
                          !grepl("_THRESH_AL",colnames(somatic))]
# Calc mean supp reads among tissues and samples by discarding from the calculation the samples where the variant is not found
somatic_supp[somatic_supp == 0] <- NA
somatic_supp$supp_reads <- apply(somatic_supp[,grepl("_ALT",colnames(somatic_supp))],1,mean,na.rm=TRUE)

# Modify names for plotting purposes
somatic_supp$somatic.anno <- gsub("_"," ",somatic_supp$somatic.anno)

# Define if common or tissue-spec
somatic_supp$type <- ifelse(grepl("Somatic common",somatic_supp$somatic.anno),"Common","Tissue-specific")

# Define if private or recurrent
somatic_supp$type2 <- ifelse(grepl("recurrent",somatic_supp$somatic.anno),"Recurrent","Private")


# Reorder
d1 <- subset(somatic_supp,type!="Common")
d2 <- subset(somatic_supp,type=="Common")
somatic_supp <- rbind(d1,d2)
somatic_supp$type <- factor(somatic_supp$type,levels = unique(somatic_supp$type))

# Dataf with variant countns
dfct <- as.data.frame(table(somatic_supp$somatic.anno))
dfct$type <- ifelse(grepl("Somatic common",dfct$Var1),"Common","Tissue-specific")
dfct$type2 <- ifelse(grepl("recurrent",dfct$Var1),"Recurrent","Private")
dfct$type <- factor(dfct$type,levels = unique(dfct$type))

# Calc median
median_df <- aggregate(somatic_supp$supp_reads,by=list(somatic_supp$somatic.anno,somatic_supp$type2),median)
names(median_df) <- c("Var1","type2","median")
median_df$type <- ifelse(grepl("Somatic common",median_df$Var1),"Common","Tissue-specific")
median_df$type2 <- ifelse(grepl("recurrent",median_df$Var1),"Recurrent","Private")
median_df$type <- factor(median_df$type,levels = unique(median_df$type))

p3=ggplot(somatic_supp,aes(x=supp_reads,y=type2,fill=type2)) +
  theme_classic() +
  geom_density_ridges(scale = .5,bandwidth=.5,alpha=.3,
                      jittered_points = F,point_size=.1) +
  coord_cartesian(xlim = c(0,10)) +
  # facet
  #facet_grid(type~.,scales = "free_y") +
  facet_grid(~type,scales = "free_y") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 10, colour = "black", face= 'bold')) +
  # x axis
  xlab("Number of supporting reads (mean)") +
  theme(axis.text.x=element_text(size = axistexty,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=axistexty+2)) +
  # y axis
  ylab("no label") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=0)) +
  # line
  #geom_vline(xintercept = c(1,2),color="grey70",linetype='dashed') +
  #geom_hline(yintercept = c(1,2)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Recurrent"="grey50","Private"="white")) +
  geom_text(data=dfct,aes(label=paste0("n=",Freq),x=8,y=type2),vjust=-5,size=3) +
  # median
  geom_text(data=median_df,aes(x=7,y=type2,label=paste0("median = ",round(median,1)," reads")),size=2.5,vjust=-2)
  
#######################################################################################################################################



#######################################################################################################################################

#### Panel D MAF

# Get reference reads
ref_reads <- str_split_fixed(somatic$REF_READS,"\\|",15)
class(ref_reads) <- "numeric"

# Get alt reads
alt_reads <- as.matrix(somatic[,grepl("_ALT",colnames(somatic)) & !grepl("_THRESH_",colnames(somatic))])

# Calc MAF (alt / (ref+alt)) for each variants in each sample
maf_mtx <- matrix(nrow = nrow(ref_reads),ncol = ncol(ref_reads))
for (i in 1:nrow(ref_reads)) {
  for (j in 1:ncol(ref_reads)) {
    ref <- ref_reads[i,j]
    alt <- alt_reads[i,j]
    # just a silly check
    if (!is.na(alt)) {if (alt>ref) {stop("ALT > REF - how is it possible?")}}
    maf <- alt / (ref+alt)
    maf_mtx[i,j] <- maf
  }
}
# Convert to df
maf_df <- as.data.frame(maf_mtx)

# Set 0 to NA
maf_df[maf_df == 0] <- NA

# Calc mean
maf_df$maf <- apply(maf_df,1,mean,na.rm=TRUE)

# Annotate
maf_df$somatic.anno <- somatic$somatic.anno

# Modify names for plotting purposes
maf_df$somatic.anno <- gsub("_"," ",maf_df$somatic.anno)

# Define if common or tissue-spec
maf_df$type <- ifelse(grepl("Somatic common",maf_df$somatic.anno),"Common","Tissue-specific")

# Define if private or recurrent
maf_df$type2 <- ifelse(grepl("recurrent",maf_df$somatic.anno),"Recurrent","Private")


# Reorder
d1 <- subset(maf_df,type!="Common")
d2 <- subset(maf_df,type=="Common")
maf_df <- rbind(d1,d2)
maf_df$type <- factor(maf_df$type,levels = unique(maf_df$type))

# Calc median
median_df2 <- aggregate(maf_df$maf,by=list(maf_df$somatic.anno,maf_df$type2),median)
names(median_df2) <- c("Var1","type2","median")
median_df2$type <- ifelse(grepl("Somatic common",median_df2$Var1),"Common","Tissue-specific")
median_df2$type2 <- ifelse(grepl("recurrent",median_df2$Var1),"Recurrent","Private")
median_df2$type <- factor(median_df2$type,levels = unique(median_df2$type))


p4=ggplot(maf_df,aes(x=maf,y=type2,fill=type2)) +
  theme_classic() +
  geom_density_ridges(scale = .5,bandwidth=.01,alpha=.3,
                      jittered_points = F,point_size=.1) +
  coord_cartesian(xlim = c(0,.25)) +
  # facet
  facet_grid(~type,scales = "free_y") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  theme(strip.text.y = element_text(size = 10, colour = "black", face= 'bold')) +
  # x axis
  xlab("MAF (mean)") +
  theme(axis.text.x=element_text(size = axistexty,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=axistexty+2)) +
  # y axis
  ylab("no label") +
  theme(axis.text.y=element_text(size = axistexty)) +
  theme(axis.title.y = element_text(size=0)) +
  # line
  #geom_vline(xintercept = c(0.05),color="grey70",linetype='dashed') +
  scale_x_continuous(breaks=c(0.05,0.15,0.25)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Recurrent"="grey50","Private"="white")) +
  geom_text(data=dfct,aes(label=paste0("n=",Freq),x=.18,y=type2),vjust=-5,size=3) +
  # median
  geom_text(data=median_df2,aes(x=0.17,y=type2,label=paste0("median MAF = ",round(median,2))),size=2.5,vjust=-2)






riga2=ggarrange(panel_B,panel_A,nrow=2,labels = c("A","B"),heights = c(.48,.52))
riga3=ggarrange(p3,p4,labels = c("C","D"),ncol=2,nrow=1)
ff=ggarrange(riga2,riga3,nrow=2,heights = c(.66,.34))
ggsave(filename='Fig2.pdf',plot=ff,width=210,height=297/2/2*3,units='mm')




