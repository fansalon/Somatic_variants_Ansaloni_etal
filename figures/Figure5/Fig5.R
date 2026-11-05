#######################################################################################################################################
### Figure - Illumina + Nanorepeat validation
#######################################################################################################################################
rm(list = ls())
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(ggh4x))
suppressMessages(library(ggridges))
suppressMessages(library(scales))
#######################################################################################################################################
# load files
smtc <- fread("../../analyses/20_Illumina_validation-sansacompvcf/ONT_somatic_validation.txt",header=T,sep='\t',data.table = F)
grml <- fread("../../analyses/20_Illumina_validation-sansacompvcf/ONT_germline_validation.txt",header=T,sep='\t',data.table = F)
# bind files
smtc$Type <- "Somatic"
grml$Type <- "Germline"
tot  <- rbind(smtc,grml)
tot$Type <- factor(tot$Type,levels = unique(tot$Type))
# remove DUP - too few
tot <- subset(tot,SVType!="DUP")
#######################################################################################################################################

#### Panel A - validation %

# count
count <- as.data.frame(table(tot[,c("Classification","Type")]))
count <- subset(count,Freq>0)

# rename variables
count$Classification <- ifelse(count$Classification=="FP","Not validated","Validated")

# validation %
totals <- aggregate(count$Freq,by=list(count$Type),sum)
names(totals) <- c("Type","tot")
totals$id <- paste(totals$Type,sep=':')
totals <- totals[,c("id","tot")]
# calc %
count$id <- paste(count$Type,sep=':')
count <- merge(count,totals,by="id",all=T,sort=F)
count$perc <- count$Freq / count$tot * 100


# plot
p1=ggplot(data=subset(count,Classification=="Validated"),aes(x=Type,y=perc,fill=Type)) +
  theme_classic() +
  geom_bar(stat='identity',color='black',alpha=.9,position = position_dodge()) +
  facet_nested(~Type,scales = "free_x",space = "free") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # x-axis
  theme(axis.text.x=element_text(angle=0, hjust=.5, size = 10)) + 
  theme(axis.title.x = element_text(size=0)) +
  # y-axis
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  ylab("% of validated variants") +
  # scale fill
  scale_fill_manual(values = c("Somatic"="white","Germline"="white")) +
  # add label
  geom_text(data=subset(count,Classification=="Validated"),
            aes(y=perc+4.5,label=paste0(round(perc,1),"%")),
            size=3,position = position_dodge(1)) +
  geom_text(data=subset(count,Classification=="Validated"),
            aes(y=perc+1.5,label=paste0(Freq,"/",tot)),size=2.5,position = position_dodge(1)) +
  theme(legend.position = 'none')

#######################################################################################################################################


#### Panel B - validation % - SR or not


# count
count2 <- as.data.frame(table(tot[,c("Classification","Type","overlaps.simplerep_ont")]))
count2 <- subset(count2,Freq>0)

# rename variables
count2$Classification <- ifelse(count2$Classification=="FP","Not validated","Validated")
count2$overlaps.simplerep_ont <- ifelse(count2$overlaps.simplerep_ont==T,
                                        "In simple repeats","Not in simple repeats")

# validation %
totals <- aggregate(count2$Freq,by=list(count2$Type,count2$overlaps.simplerep_ont),sum)
names(totals) <- c("Type","overlaps.simplerep_ont","tot")
totals$id <- paste(totals$Type,totals$overlaps.simplerep_ont,sep=':')
totals <- totals[,c("id","tot")]
# calc %
count2$id <- paste(count2$Type,count2$overlaps.simplerep_ont,sep=':')
count2 <- merge(count2,totals,by="id",all=T,sort=F)
count2$perc <- count2$Freq / count2$tot * 100


# plot
p2=ggplot(data=subset(count2,Classification=="Validated"),aes(x=Type,y=perc,fill=overlaps.simplerep_ont)) +
  theme_classic() +
  geom_bar(stat='identity',color='black',alpha=.7,position = position_dodge()) +
  facet_nested(~Type,scales = "free_x",space = "free") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # x-axis
  theme(axis.text.x=element_text(angle=0, hjust=.5, size = 10)) + 
  theme(axis.title.x = element_text(size=0)) +
  # y-axis
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  ylab("% of validated variants") +
  # scale fill
  scale_fill_manual(values = c("In simple repeats"="darkgoldenrod1","Not in simple repeats"="grey75")) +
  # add label
  geom_text(data=subset(count2,Classification=="Validated"),
            aes(y=perc+6,label=paste0(round(perc,1),"%")),
            size=3,position = position_dodge(1)) +
  geom_text(data=subset(count2,Classification=="Validated"),
            aes(y=perc+2,label=paste0(Freq,"/",tot)),size=2.2,position = position_dodge(1)) +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(title=""))




#######################################################################################################################################


#### Panel C - validation % - per svtype + SR or not


# count
count1 <- as.data.frame(table(tot[,c("Classification","SVType","Type","overlaps.simplerep_ont")]))
count1 <- subset(count1,Freq>0)

# rename variables
count1$Classification <- ifelse(count1$Classification=="FP","Not validated","Validated")
count1$overlaps.simplerep_ont <- ifelse(count1$overlaps.simplerep_ont==T,
                                        "In simple repeats","Not in simple repeats")

# validation %
totals <- aggregate(count1$Freq,by=list(count1$SVType,count1$Type,count1$overlaps.simplerep_ont),sum)
names(totals) <- c("SVType","Type","overlaps.simplerep_ont","tot")
totals$id <- paste(totals$SVType,totals$Type,totals$overlaps.simplerep_ont,sep=':')
totals <- totals[,c("id","tot")]
# calc %
count1$id <- paste(count1$SVType,count1$Type,count1$overlaps.simplerep_ont,sep=':')
count1 <- merge(count1,totals,by="id",all=T,sort=F)
count1$perc <- count1$Freq / count1$tot * 100

# validation %, per svtype
lbl <- aggregate(count1$Freq,by=list(count1$Classification,count1$SVType,count1$Type),sum)
names(lbl) <- c("Classification","SVType","Type","Freq")
lbl$tot <- aggregate(count1$tot,by=list(count1$Classification,count1$SVType,count1$Type),sum)[,"x"]
lbl$perc <- lbl$Freq / lbl$tot *100

# plot
p3=ggplot(data=subset(count1,Classification=="Validated"),aes(x=SVType,y=perc)) +
  theme_classic() +
  geom_bar(stat='identity',color='black',alpha=.7,position = position_dodge(),aes(fill=overlaps.simplerep_ont)) +
  facet_nested(~Type,scales = "free_x",space = "free") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # x-axis
  theme(axis.text.x=element_text(angle=0, hjust=.5, size = 10)) + 
  theme(axis.title.x = element_text(size=0)) +
  # y-axis
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  ylab("% of validated variants") +
  # scale fill
  scale_fill_manual(values = c("In simple repeats"="darkgoldenrod1","Not in simple repeats"="grey75")) +
  # add label
  geom_text(data=subset(count1,Classification=="Validated"),
            aes(y=perc+6,label=paste0(round(perc,1),"%"),fill=overlaps.simplerep_ont),
            size=3,position = position_dodge(1)) +
  geom_text(data=subset(count1,Classification=="Validated"),
            aes(y=perc+2,label=paste0(Freq,"/",tot),fill=overlaps.simplerep_ont),size=2.2,position = position_dodge(1)) +
  geom_text(data=subset(lbl,Classification=="Validated"),
            aes(y=53,label=paste0(round(perc,1),"%")),size=3,position = position_dodge(1)) +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(title=""))

#######################################################################################################################################

#######################################################################################################################################


### Panels D and E

tot$Classification <- ifelse(tot$Classification=="TP","Validated",
                             ifelse(tot$Classification=="FP","Not validated","ERR"))

p4=ggplot(tot,aes(y=Type,x=MAF_ont,color=Classification,fill=Classification)) +
  theme_classic() +
  geom_density_ridges(scale = .5,bandwidth=.01,alpha=.2,
                      jittered_points = F,point_size=.1) +
  # facet
  #facet_nested(~Classification) +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # x axis
  xlab("MAF - ONT variants") +
  theme(axis.text.x=element_text(size = 10,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=11)) +
  # y axis
  ylab("no label") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(legend.title = element_blank()) +
  # change color
  scale_color_manual(values = c("Validated"="#00BA38","Not validated"="red")) +
  scale_fill_manual(values = c("Validated"="#00BA38","Not validated"="red"))
#  theme(legend.position = 'none')
p5=ggplot(tot,aes(y=Type,x=ALT_reads_ont,color=Classification,fill=Classification)) +
  theme_classic() +
  geom_density_ridges(scale = .5,bandwidth=.08,alpha=.2,
                      jittered_points = F,point_size=.1) +
  # facet
  # facet_nested(~Classification) +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # x axis
  xlab("Supporting reads - ONT variants") +
  theme(axis.text.x=element_text(size = 10,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=11)) +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  # y axis
  ylab("no label") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(legend.title = element_blank()) +
  # change color
  scale_color_manual(values = c("Validated"="#00BA38","Not validated"="red")) +
  scale_fill_manual(values = c("Validated"="#00BA38","Not validated"="red"))
#  theme(legend.position = 'none')

#######################################################################################################################################

#######################################################################################################################################


### Panels F - variants MAF>0.99


germline_maf9 <- subset(tot,MAF_ont>0.99)
counts <- as.data.frame(table(germline_maf9[,c("Classification")]))

p6=ggplot(counts,aes(x="Variants\nMAF > 0.99",y=Freq)) +
  theme_classic() +
  geom_bar(stat='identity',aes(fill=Var1),color='black',position = position_dodge(.9),alpha=.7) +
  # x axis
  xlab("Supporting reads - ONT variants") +
  theme(axis.text.x=element_text(size = 10,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("No. of variants MAF > 0.99") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(legend.title = element_blank()) +
  geom_text(aes(label=Freq,fill=Var1,y=Freq+100),size=3,position = position_dodge(.9)) +
  # change fill color
  scale_fill_manual(values = c("Validated"="#00BA38","Not validated"="red"))
#######################################################################################################################################

#######################################################################################################################################


### Panels H - nanorepeat results

# result dir
resdir <- "../../analyses/21_nanorepeat_validation/output/"
sd_interval <- 3


### Load results
tissues <- c("CTX","HIP","KID")
mouse   <- c(2,3,4,5,6)
i <- 0
for (tx in tissues) {
  for (mm in mouse) {  
    i <- i + 1
    
    ### summary
    df1 <- fread(paste0(resdir,"/summary_WT_",mm,"_",tx,".txt"),header=T,sep='\t',data.table = F)
    df1$Validated <- ifelse(df1$Validated==FALSE & df1$no_GT=="1","FALSE-no ALT genotype",
                            ifelse(df1$Validated==FALSE & df1$no_GT=="2","FALSE-ALT genotype supports INS",df1$Validated))
    
    
    ### validation results
    vr1 <- fread(paste0(resdir,"/validation_res_WT_",mm,"_",tx,".txt"),header=T,sep='\t',data.table = F)
    
    # define interval SVLEN +/-3sd and the same for nanorepeat DEL
    vr1$SVLEN_start <- vr1$SVLEN - vr1$SVLEN_sd*sd_interval
    vr1$SVLEN_end <- vr1$SVLEN + vr1$SVLEN_sd*sd_interval
    vr1$DEL_len_start <- abs(vr1$DEL_len) - vr1$DEL_len_sd*sd_interval
    vr1$DEL_len_end <- abs(vr1$DEL_len) + vr1$DEL_len_sd*sd_interval
    
    # define if SVLEN +/- 3 sd overlaps nanorepeat DEL +/- 3 sd
    vr1$val <- ifelse(vr1$SVLEN_start <= vr1$DEL_len_end & vr1$DEL_len_start <= vr1$SVLEN_end,T,F)
    
    # get number of negative events
    neg <- table(vr1$val)["FALSE"]
    
    ### Update the summary
    df1$Count <- ifelse(df1$Validated==TRUE,df1$Count-neg,df1$Count)
    new_line <- data.frame(no_GT=2,Count=neg,Supported_SVTYPE="DEL",
                           Validated="FALSE-ALT genotype supports DEL of different length",
                           Sample=unique(df1$Sample))
    df2 <- rbind(df1,new_line)
    rownames(df2) <- NULL
    
    # bind
    if (i==1) {
      summary_tot <- df2
    } else {
      summary_tot <- rbind(summary_tot,df2)
    }
  }
}


### Plot results

# discard "no" and "multiple" genotypes - this have not been analysed
summary_tot <- subset(summary_tot,no_GT!="0" & no_GT!=">2")

# calc %
tot <- aggregate(summary_tot$Count,by=list(summary_tot$Sample),sum)
names(tot) <- c("Sample","tot")
summary_tot <- merge(summary_tot,tot,by="Sample",all=T,sort=F)
summary_tot$Perc <- summary_tot$Count / summary_tot$tot *100

# define tissue and mouse number
summary_tot$Tissue <- str_split_fixed(summary_tot$Sample,"_",3)[,3]
summary_tot$mouse <- paste0("Mouse ",as.numeric(str_split_fixed(summary_tot$Sample,"_",3)[,2])-1)

p7=ggplot(summary_tot,aes(x=Tissue,y=Count,fill=Validated)) +
  theme_classic() +
  # bars
  geom_bar(stat='identity',color='black',alpha=.8) +
  # text
  geom_text(data=subset(summary_tot,Validated==TRUE),
            size =3, position = position_stack(vjust = .9),aes(label=paste0(round(Perc),"%"))) +
  # scale colors
  scale_fill_manual(values = c("FALSE-no ALT genotype"="red",
                               "FALSE-ALT genotype supports INS"="firebrick",
                               "FALSE-ALT genotype supports DEL of different length"="pink",
                               "TRUE"="#00BA38")) +
  # facet
  facet_grid(~mouse,scales = "free") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face= 'bold')) +
  # x axis
  xlab("none") +
  theme(axis.text.x=element_text(size = 12,angle = 0, hjust = .5)) +
  theme(axis.title.x = element_text(size=0)) +
  # y axis
  ylab("No. of somatic DELs\nin simple repeats (Sniffles2)") +
  theme(axis.text.y=element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) +
  # legend
  guides(fill=guide_legend(nrow=2,byrow=TRUE,title="Validated by\nNanoRepeat?")) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(.1, 'cm'),
        legend.spacing.y = unit(.001, 'cm'),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) + 
  # y limit
  coord_cartesian(ylim = c(0,6500))
#######################################################################################################################################

#######################################################################################################################################



riga1_col1=ggarrange(p1,nrow=2,labels = "A",heights = c(.88,.12))
riga1_col23=ggarrange(p2,p3,ncol=2,common.legend = T,legend = 'bottom',widths = c(.33,.66),labels = c("B","C"))
riga1=ggarrange(riga1_col1,riga1_col23,ncol=2,widths = c(.21,.79))
riga2=ggarrange(p4,p5,p6,ncol = 3,labels = c("D","E","F"),common.legend = T,legend = 'bottom',
                  widths = c(.39,.39,.22))
ff=ggarrange(riga1,riga2,ggarrange(p7,labels = "G"),nrow = 3)
ggsave(filename='Fig5.pdf',plot=ff,width=230,height=250,units='mm')






