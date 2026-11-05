library(data.table)

ff=fread("input/tldr_res.txt",data.table=F,sep='\t',header=T)


# Select only MEI: PASS, with TSD and not already known NonRef insertions
filtered <- subset(ff,Filter=="PASS" & !is.na(TSD) & is.na(NonRef) )

write.table(filtered,"tldr_res_filtered.txt",row.names=F,sep='\t',quote=F)
