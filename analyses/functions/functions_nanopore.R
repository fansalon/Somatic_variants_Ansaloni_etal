### ---------------------------------------------------------------------------
### g:Profiler version pinning
###
### gost() queries a LIVE server whose annotation database changes over time.
### Re-running against the default URL silently changes every p-value.
### Archived instances used in this paper:
###   e111_eg58_p18 - analyses 05, 08 (run Dec 2024), GO:BP releases/2024-01-17
###   e112_eg59_p19 - analysis 15     (run Apr 2025)
### ---------------------------------------------------------------------------

gp_archive_url <- function(gp_version) {
  paste0("https://biit.cs.ut.ee/gprofiler_archive3/", gp_version)
}

# hash for the result cache; uses digest if installed, base R otherwise
.gp_key <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(x))
  f <- tempfile(); on.exit(unlink(f))
  saveRDS(x, f, compress = FALSE, version = 2)
  unname(tools::md5sum(f))
}
doGprof2 <- function(query_gene, bg_genes, fdr, GO_name, n_significant_genes,
                     gp_version,                        # REQUIRED, e.g. "e111_eg58_p18"
                     cache_dir = "gprofiler_cache",
                     max_tries = 5) {

  if (missing(gp_version) || !nzchar(gp_version))
    stop("gp_version must be given explicitly (e.g. 'e111_eg58_p18')")
  gprofiler2::set_base_url(gp_archive_url(gp_version))

  # cache lookup 
  key <- .gp_key(list(gp_version = gp_version,
                      query      = sort(unique(query_gene)),
                      background = sort(unique(bg_genes)),
                      fdr        = fdr,
                      source     = GO_name))
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_file <- file.path(cache_dir, paste0(gp_version, "_", key, ".Rds"))

  if (file.exists(cache_file)) {
    gostres <- readRDS(cache_file)
    cat("\t\t[cached]\n")
  } else {
    # query, retrying on transient failures (the archives return 504)
    gostres <- NULL
    for (i in seq_len(max_tries)) {
      attempt <- tryCatch(
        gost(query             = query_gene,
             custom_bg         = bg_genes,
             organism          = "mmusculus",
             significant       = TRUE,
             user_threshold    = fdr,
             correction_method = "g_SCS",
             sources           = GO_name,
             exclude_iea       = TRUE,
             domain_scope      = "custom_annotated",   # FONDAMENTALE
             evcodes           = TRUE),
        error = function(e) e)
      if (!inherits(attempt, "error")) { gostres <- attempt; break }
      if (i == max_tries) stop(attempt)
      wait <- 120 * i
      cat("\t\tg:Profiler failed (", conditionMessage(attempt),
          ") - retry ", i, "/", max_tries - 1, " in ", wait, "s\n", sep = "")
      Sys.sleep(wait)
    }
    # confirm we really reached the pinned instance
    if (!is.null(gostres)) {
      got <- gostres$meta$version
      if (!startsWith(got, gp_version))
        stop("Requested g:Profiler ", gp_version, " but the server answered as ", got)
    }
    saveRDS(gostres, cache_file)      # NULL (= no significant terms) is cached too
  }

  if (is.null(gostres)) {
    res_filt_sub <- data.frame()
  } else {
    res <- gostres$result
    res_filt <- subset(res, intersection_size >= n_significant_genes)
    res_filt_sub <- res_filt[, c("term_id","term_name","intersection_size","query_size",
                                 "term_size","effective_domain_size","p_value","intersection")]
    names(res_filt_sub)[ncol(res_filt_sub)]     <- "Gene.name"
    names(res_filt_sub)[ncol(res_filt_sub) - 1] <- "p.adj"
  }
  return(res_filt_sub)
}
getNumbReads <- function(df,cols,i) {
  df <- df[i,]
  # get the number of ALT reads for the samples 0/0
  alt_0.0 <- df[, gsub("_GT","_ALT",cols[ df[,cols]=="0/0" ]) ] 
  # if not a dataf, it means only one 0/0 sample was found
  if ( is.data.frame(alt_0.0) ) {
    counter <- rowSums ( alt_0.0 >= sup_reads )
  } else {
    if (alt_0.0 >= sup_reads) {
      counter <- 1
    } else {
      counter <- 0
    }
  }
  return(counter)
}
extractFeatureFromTxDb <- function(TxDb,upstreamTSS,downstreamTSS,downstreamTTS) {
  #### Retrieve genomic features from txdb
  #print("Retrieving coordinates of the genomic features")
  ## Transcripts
  transcripts <- transcripts(TxDb,columns=c("tx_name","gene_id"),use.names=T)
  ## Promoters
  promoters <- promoters(TxDb, upstream=upstreamTSS, downstream=downstreamTSS,columns=c("tx_name","gene_id"),use.names=T)
  ## Exons
  exons <- exons(TxDb,use.names=T,columns=c("tx_name","gene_id","exon_name"))
  # expand exons (tx_name and gene_id is actually a list of transcripts sharing the same exons)
  exons <- S4Vectors::expand(exons, c("tx_name","gene_id"))
  mcols(exons)$tx_name <- exons$tx_name
  ## 5 and 3' utrs
  # utr5 and utr3 are genomic ranges list, convert them to GR with unlist and add transcript name as metadata with mcols
  utr5 <- unlist(fiveUTRsByTranscript(TxDb,use.names=T))
  mcols(utr5)$tx_name <- names(utr5)
  utr3 <- unlist(threeUTRsByTranscript(TxDb,use.names=T))
  mcols(utr3)$tx_name <- names(utr3)
  ## Introns
  introns <- unlist(intronsByTranscript (TxDb,use.names=TRUE))
  mcols(introns)$tx_name <- names(introns)
  ## Downstream (transcript TTS + downstreamTTS)
  downstream <- flank(transcripts, downstreamTTS, start=FALSE)
  ### combine everything in a GR list
  grl <- GRangesList(Promoter=promoters,Exon=exons,utr5=utr5,utr3=utr3,Intron=introns,Downstream=downstream)
  return(grl)
}
getSampleSV <- function(df,reads_thresh) {
  i <- 0
  for (j in 1:ncol(df)) {
    name <- strsplit(colnames(df)[j],"_")[[1]]
    if (name[1]=="WT" & (name[2]=="1"|name[2]=="2"|name[2]=="3"|name[2]=="4"|name[2]=="5"|name[2]=="6") & 
        name[4]=="ALT") {
      i <- i + 1
      # Get SV supported by at least reads_thresh in this sample
      new <- subset(df,df[,j]>=reads_thresh)
      # Get sample name
      sample <- paste(name[1],name[2],name[3],sep="_")
      new$Sample <- sample
      # Get GENOTYPE.QUAL
      new$GENOTYPE.QUAL <- as.numeric(str_split_fixed(new$GENOTYPE.QUAL_VEC,"\\|",18)[,i])
      # Get REF reads
      new$REF_READS <- as.numeric(str_split_fixed(new$REF_READS,"\\|",18)[,i])
      # Get genotype
      new$Genotype <- new[,grepl(paste0(sample,"_GT"),colnames(new))]
      # Get ALT reads
      new$ALT_READS <- new[,grepl(paste0(sample,"_ALT"),colnames(new))]
      # Set sample
      new$Sample <- sample
      # Get columns of interest
      new.f <- new[,c("chrom","start","end","FILTER","PRECISION","SVTYPE","SVLEN","STRAND","breakpoint2","GENOTYPE.QUAL","REF_READS","ALT_READS","Genotype","Sample")]
      # Set rep number
      new.f$Rep <- paste0("Rep_",name[2])
      # Set tissue
      new.f$Tissue <- name[3]
      # bind
      if (i==1) {
        final <- new.f
      } else {
        final <- rbind(final,new.f)
      }
    }
  }
  return(final)
}
annotateGenomicFeatures <- function(df,priorities,txdbGR,quiet,is.shuffle,txdb_raw) {
  
  # create SV id and convert file to GR
  df$id <- paste0(df$chrom,"-",df$start,"-",df$end,"-",df$Sample)
  
  # get genomic features
  feature_Grl <- txdbGR
  dfgr2 <- annotateGR(df=df,TxDb=feature_Grl,quiet,prior=priorities,is.shuffle,txdb_raw)
  
  return(dfgr2)
}
annotateGR <- function(df,TxDb,quiet,prior,is.shuffle,txdb_raw) {
  ### Define priorities
  if (length(prior) != 6) {
    stop("the vector of priorities does not contain 6 features")
  }
  if (! quiet) {
    print(paste0("Assigning  priorities as follow: ",prior[1],", ",prior[2],", ",prior[3],", ",prior[4],", ",
                 prior[5],", ",prior[6]))
  }
  
  ### Convert to granges
  gr_bed <- toGRanges(df)
  
  ### Overlap gr_bed with feature_Grl - return is dataframe (chr,start,end,id,txname,annotation)
  overlaps <- grOverlapTxDb(GR_bed=gr_bed,grl_txdb=TxDb,c_prior=prior,quiet=quiet,is.shuffle=is.shuffle,txdb=txdb_raw)
  
  ### Add the annotation to the original file
  if (is.shuffle) {
    overlaps <- overlaps[,c("id","annotation")]
  } else {
    overlaps <- overlaps[,c("id","tx_name2","tx_coord","distanceTSS","distanceTx","annotation")]
    colnames(overlaps)[colnames(overlaps)=="tx_name2"] <- "tx_name"
  }
  final <- merge(df,overlaps,by="id",sort=F)
  
  
  ### Check input and output have the same length
  if (nrow(final)!=nrow(df)) {
    stop("input and output do not have the same length")
  }
  return(final)
}
grOverlapTxDb <- function(GR_bed,grl_txdb,c_prior,quiet,is.shuffle,txdb) {
  
  nro1 <- length(GR_bed)
  counter <- 0
  ### Start to loop among the priorities (promoter, utr, etc, etc)
  for (prrr in c_prior) {
    counter <- counter + 1
    if (! quiet ) { print(paste0("Overlapping with ",prrr)) }
    # select the coords from the txdb (promoter, utr, etc, etc)
    txdb_gr <- grl_txdb[[prrr]]
    # Overlap the 2 GR files (e.g. recomb with promoter coords)
    ovrlp <- suppressWarnings(findOverlaps(query=GR_bed, subject=txdb_gr))
    # Extract row indices for query (recomb) and subject (promoter)
    quer <- as.data.frame(GR_bed[queryHits(ovrlp)],row.names = NULL)
    sub <- as.data.frame(txdb_gr[subjectHits(ovrlp)],row.names = NULL)
    # rename in order not to have duplicated names
    names(sub) <- paste0(names(sub),"2")
    # combine the 2 dataframes
    merged <- unique(cbind(quer,sub)[,c("seqnames","start","end","id","start2","end2","strand2","tx_name2","width2")])
    
    if (nrow(merged)) {
      
      # if analysing real and not shuffle SV, report also associated transcript
      if (!is.shuffle) {
        
        # calc distance from TSS (if feature is promoter)
        if (prrr=="Promoter") {
          
          # identify TSS (which is the middle point in between the promoter interval)
          prom_leng <- unique(merged$width2)/2
          merged$TSS <- ifelse(merged$strand2=="+",merged$start2+prom_leng,merged$end2-prom_leng)
          
          # calc distance to TSS
          merged$distanceTSS_1 <- abs(merged$start-merged$TSS) # distance from TSS to start of SV
          merged$distanceTSS_2 <- abs(merged$end-merged$TSS)  # distance from TSS to end of SV
          merged$distanceTSS <- apply(merged[,c("distanceTSS_1","distanceTSS_2")],1,min) # real distance is the minimum of the two
          # if SV contains TSS, distance is 0
          merged$distanceTSS <- ifelse(merged$start<=merged$TSS & merged$end>=merged$TSS,0,merged$distanceTSS)
          merged <- merged[,c("seqnames","start","end","id","tx_name2","distanceTSS")]
          
          # handle redundant overlaps (closest to TSS, if more then one closest, get one random?)
          merged <- as.data.frame(merged %>% 
                                    group_by(id) %>% 
                                    slice_min(distanceTSS, n = 1,with_ties=F) )
          
        } else {
          
          # distance to TSS is set to NA
          merged$distanceTSS <- NA
          
          # calc length of overlap
          merged$overlap <- 
            # means that the feature is fully contained in the SV
            ifelse( merged$start<=merged$start2 & merged$end>=merged$end2,merged$end2-merged$start2,
                    # means that the SV is fully contained in the feature
                    ifelse(merged$start>=merged$start2 & merged$end<=merged$end2,merged$end-merged$start,
                           # means that the SV is upstream to the feature
                           ifelse(merged$start<merged$start2 & merged$end<merged$end2,merged$end-merged$start2,
                                  # means that the SV is downstream to the feature
                                  ifelse(merged$start>merged$start2 & merged$end>merged$end2,merged$end2-merged$start,NA))))
          merged$overlap <- ifelse(merged$overlap==0,1,merged$overlap) # it can be distance is 0 if some SV have weirdly len=0 or SV and feature are adjacent
          
          # # test start
          # unique(merged[merged$start<=merged$start2 & merged$end>=merged$end2 ,c("seqnames","start","end","start2","end2","overlap")])
          # unique(merged[merged$start>=merged$start2 & merged$end<=merged$end2 ,c("seqnames","start","end","start2","end2","overlap")])
          # unique(merged[merged$start<merged$start2 & merged$end<merged$end2 ,c("seqnames","start","end","start2","end2","overlap")])
          # unique(merged[merged$start>merged$start2 & merged$end>merged$end2 ,c("seqnames","start","end","start2","end2","overlap")])
          # # test end
          
          # handle redundant overlaps (select the max overlap, or random)
          merged <- as.data.frame(merged %>% 
                                    group_by(id) %>% 
                                    slice_max(overlap, n = 1,with_ties=F) )
          merged <- merged[,c("seqnames","start","end","id","tx_name2","distanceTSS")]
          
        }
      } else {
        merged$tx_name2 <- NA
        merged$distanceTSS <- NA
        merged <- unique(merged[,c("seqnames","start","end","id","tx_name2","distanceTSS")])
      }
    }
    
    
    #merged$gene_id <- unlist(merged$gene_id)
    # Add genomic feature name - if not empty
    if (nrow(merged)>0) {
      merged$distanceTx <- 0 # to make it compatible with intergenic SVs
      merged$annotation <- prrr
      ##############################################################################################################################
      ####### # select non redundant terms (when redundant ids, the first is selected)
      #######merged_nred <- merged %>% distinct(id, .keep_all = T)
      # discard from the GR_bed the intervals that were positively associated to the analysed feature
      ##############################################################################################################################
      positive <- merged$id
      GR_bed <- GR_bed[!GR_bed$id %in% positive]
      # write to final file
      if (counter == 1) {
        final <- merged
      } else {
        final <- rbind(merged,final)
      }
    } else {
      if (counter == 1) {
        final <- data.frame()
      }
    }
  }
  
  # if not empty
  if (length(GR_bed)>0) {
    # if not shuffle
    if (!is.shuffle) {
      ### Finally, the un-annotated ones are annotated as Intergenic
      # Get closest transcript
      dist_interg <- distanceToNearest(GR_bed,transcripts(txdb),select="all")
      # Extract row indices for query (SV) and subject (tx)
      quer <- as.data.frame(GR_bed[queryHits(dist_interg)],row.names = NULL)
      sub <- as.data.frame(txdb_gr[subjectHits(dist_interg)],row.names = NULL)
      # rename in order not to have duplicated names
      names(sub) <- paste0(names(sub),"2")
      # combine the 2 dataframes
      merged <- (cbind(quer,sub)[,c("seqnames","start","end","id","tx_name2")])
      merged$distanceTSS <- NA
      # add distance tot tx
      merged$distanceTx <- as.data.frame(dist_interg)[,"distance"]
      # handle redundant overlaps (select closest tx, or random)
      merged <- as.data.frame(merged %>% 
                                group_by(id) %>% 
                                slice_min(distanceTx, n = 1,with_ties=F) )
      # check non-redundant selection worked
      if (nrow(merged) != length(unique(merged$id))) { stop("found redundant intergenic SVs")}
      # annotate
      merged$annotation <- "Distal Intergenic"
      # bind with final
      final <- rbind(merged,final)
    } else {
      merged <- as.data.frame(GR_bed)[,c("seqnames","start","end","id")]
      merged$tx_name2 <- NA
      merged$distanceTSS <- NA
      merged$distanceTx <- NA
      merged$annotation <- "Distal Intergenic"
      final <- rbind(merged,final)
    }
  }
  # just check input and output have the same number of rows
  if (nrow(final) != nro1) { stop("df before and after annotation does not have the same number of lines")}
  
  
  
  ### Add coord tx
  if (!is.shuffle) {
    all_tx <- as.data.frame(transcripts(txdb))
    all_tx$tx_coord <- paste0(all_tx$seqnames,":",all_tx$start,"-",all_tx$end)
    all_tx <- all_tx[,c("tx_coord","tx_name")]
    # check all tx of interest are in all_tx
    if (!table(final$tx_name2 %in% all_tx$tx_name)["TRUE"]==nrow(final)) { stop("Some transcripts in final are missing from all_tx") }
    final2 <- merge(final,all_tx,by.x="tx_name2",by.y="tx_name",sort=F,all.x=T)
    # reorder cols
    final2 <- final2[,c( colnames(final),"tx_coord") ]
    
    # just check input and output have the same number of rows
    if (nrow(final2) != nro1) { stop("df before and after annotation does not have the same number of lines")}
  } else {
    final2 <- final
  }
  
  return(final2)
}
annotateGenomicFeaturesBND <- function(df,priorities,txdbGR,quiet,is.shuffle,txdb_raw) {
  
  # create SV id and convert file to GR
  df$id <- paste0(df$chrom,"-",df$start,"-",df$end,"-",df$sample_comparison)
  
  # get genomic features
  feature_Grl <- txdbGR
  dfgr2 <- annotateGR_BND(df=df,TxDb=feature_Grl,quiet,prior=priorities,is.shuffle,txdb_raw)
  
  return(dfgr2)
}
annotateGR_BND <- function(df,TxDb,quiet,prior,is.shuffle,txdb_raw) {
  ### Define priorities
  if (length(prior) != 6) {
    stop("the vector of priorities does not contain 6 features")
  }
  if (! quiet) {
    print(paste0("Assigning  priorities as follow: ",prior[1],", ",prior[2],", ",prior[3],", ",prior[4],", ",
                 prior[5],", ",prior[6]))
  }
  
  ######## Breakpoint 1
  
  ### Convert to granges
  gr_bed_1 <- toGRanges(df)
  ### Overlap gr_bed with feature_Grl - return is dataframe (chr,start,end,id,txname,annotation)
  overlaps_1 <- grOverlapTxDb(GR_bed=gr_bed_1,grl_txdb=TxDb,c_prior=prior,quiet=quiet,is.shuffle=is.shuffle,txdb=txdb_raw)
  
  
  ######## Breakpoint 2
  
  ### Convert to granges
  df2 <- df
  df2$chrom <- str_split_fixed(df2$breakpoint2,":",2)[,1]
  df2$start <- as.numeric(str_split_fixed(df2$breakpoint2,":",2)[,2])
  df2$end <- df2$start + 1
  
  gr_bed_2 <- toGRanges(df2)
  ### Overlap gr_bed with feature_Grl - return is dataframe (chr,start,end,id,txname,annotation)
  overlaps_2 <- grOverlapTxDb(GR_bed=gr_bed_2,grl_txdb=TxDb,c_prior=prior,quiet=quiet,is.shuffle=is.shuffle,txdb=txdb_raw)
  
  
  ######## Merge Breakpoint 1 and 2
  if ( table(overlaps_1$id %in% overlaps_2$id)["TRUE"]!=nrow(overlaps_1)){ stop() }
  overlaps <- merge(overlaps_1,overlaps_2,by="id")
  overlaps$annotation <- ifelse(overlaps$annotation.x >= overlaps$annotation.y,
                                paste0(overlaps$annotation.x,"/",overlaps$annotation.y),
                                paste0(overlaps$annotation.y,"/",overlaps$annotation.x))
  
  
  ### Add the annotation to the original file
  if (is.shuffle) {
    overlaps <- overlaps[,c("id","annotation")]
  } else {
    overlaps_1 <- overlaps_1[,c("id","tx_name2","tx_coord","distanceTSS","distanceTx","annotation")]
    overlaps_2 <- overlaps_2[,c("id","tx_name2","tx_coord","distanceTSS","distanceTx","annotation")]
    # unique annotation formatting
    overlaps$tx_name2 <- ifelse(overlaps$annotation.x >= overlaps$annotation.y,
                                paste0(overlaps$tx_name2.x,"/",overlaps$tx_name2.y),
                                paste0(overlaps$tx_name2.y,"/",overlaps$tx_name2.x))
    
    overlaps$tx_coord <- ifelse(overlaps$annotation.x >= overlaps$annotation.y,
                                paste0(overlaps$tx_coord.x,"/",overlaps$tx_coord.y),
                                paste0(overlaps$tx_coord.y,"/",overlaps$tx_coord.x))
    
    overlaps$distanceTSS <- ifelse(overlaps$annotation.x >= overlaps$annotation.y,
                                   paste0(overlaps$distanceTSS.x,"/",overlaps$distanceTSS.y),
                                   paste0(overlaps$distanceTSS.y,"/",overlaps$distanceTSS.x))
    
    overlaps$distanceTx <- ifelse(overlaps$annotation.x >= overlaps$annotation.y,
                                  paste0(overlaps$distanceTx.x,"/",overlaps$distanceTx.y),
                                  paste0(overlaps$distanceTx.y,"/",overlaps$distanceTx.x))
    overlaps <- overlaps[,c("id","annotation","tx_name2","tx_coord","distanceTSS","distanceTx")]
    colnames(overlaps)[colnames(overlaps)=="tx_name2"] <- "tx_name"
  }
  final <- merge(df,overlaps,by="id",sort=F)
  
  
  ### Check input and output have the same length
  if (nrow(final)!=nrow(df)) {
    stop("input and output do not have the same length")
  }
  return(final)
}
countTEoverlap2 <- function(path_dataf,dataf_melt,rmsk,bedtools,i) {
  
  # intersect SV with TE coords
  cmd <- paste0(bedtools," intersect -wao -f .5 -a ",path_dataf," -b ",rmsk," > tmp.intersect_bootstr",i)
  #cat("Executing",cmd,"\n")
  system(cmd)
  # read results
  intersect.orig <- fread(paste0("tmp.intersect_bootstr",i),header=F,sep='\t',data.table = F)
  system(paste0("rm tmp.intersect_bootstr",i))
  
  # assign column names
  col1 <- colnames(fread(path_dataf,header=T,sep='\t'))
  col2 <- c("chr_TE","start_TE","end_TE","TE1","TE2","strand_TE","intersect")
  header1 <- c(col1,col2)
  if (length(header1)!=ncol(intersect.orig)) {stop("error in the column formatting")}
  colnames(intersect.orig) <- header1
  
  # select columns of interest
  intersect <- unique(intersect.orig[,c("ID","TE1","TE2","intersect")])
  intersect$TE <- paste0(intersect$TE1,":",intersect$TE2)
  names(intersect) <- c("ID","discard1","discard2","overlap","TE")
  intersect$discard1 <- NULL
  intersect$discard2 <- NULL
  # assign -1 to SV not overlapping TE
  intersect$TE <- ifelse(intersect$TE==".:-1","No TE",intersect$TE)
  # remove redundancies (one SV can overlap more than one TE)
  selected <- as.data.frame(intersect %>% 
                              group_by(ID) %>% 
                              slice_max(overlap, n = 1,with_ties=F) )
  if (nrow(selected)!=length(unique(intersect$ID))) { stop("Redundant SVs found") }
  # Add info to melted file
  selected <- selected[,c("ID","TE")]
  uno <- dataf_melt
  uno <- merge(uno,selected,by="ID",all.x=T,sort=F)
  # sanity check
  if (nrow(uno[is.na(uno$TE),])>0) { stop("Some SVs have not been associated to any TE") }
  ### count
  uno$TE2 <- ifelse(uno$TE=="No TE","No TE",str_split_fixed(uno$TE,":",3)[,3])
  res2 <- as.data.frame(table(uno[,c("TE2","Sample")]))
  
  out_ls <- list(res2,uno)
  
  return(out_ls)
}
countSRoverlap <- function(path_dataf,rmsk,bedtools,i) {
  
  # intersect SV with TE coords
  cmd <- paste0(bedtools," intersect -wo -f .2 -a ",path_dataf," -b ",rmsk," > tmp.intersect_bootstr",i)
  #cat("Executing",cmd,"\n")
  system(cmd)
  # read results
  intersect.orig <- fread(paste0("tmp.intersect_bootstr",i),header=F,sep='\t',data.table = F)
  system(paste0("rm tmp.intersect_bootstr",i))
  
  # assign column names
  #col1 <- colnames(fread(path_dataf,header=T,sep='\t'))
  col1 <- c("chrom","start","end","ID")
  col2 <- c("chr_SR","start_SR","end_SR","SR","SR2","strand_SR","intersect")
  header1 <- c(col1,col2)
  if (length(header1)!=ncol(intersect.orig)) {stop("error in the column formatting")}
  colnames(intersect.orig) <- header1
  
  # get len and ID
  intersect.orig$len <- intersect.orig$end_SR-intersect.orig$start_SR
  intersect.orig$SR_ID <- paste(intersect.orig$chr_SR,intersect.orig$start_SR,intersect.orig$end_SR,sep='-')
  
  # select columns of interest
  intersect <- unique(intersect.orig[,c("ID","SR","intersect","len","SR_ID")])
  # remove redundancies (one SV can overlap more than one TE)
  selected <- as.data.frame(intersect %>% 
                              group_by(ID) %>% 
                              slice_max(order_by = tibble(intersect,len), n = 1,with_ties=F) ) # first select the one with the longest intersect, if intersection length is the same, get the variant overlapping the longest SR
  if (nrow(selected)!=length(unique(intersect$ID))) { stop("Redundant SVs found") }
  # Add info to melted file
  selected <- selected[,c("ID","SR","intersect","SR_ID")]
  
  return(selected)
}
getsample <- function(df1) {
  
  # get only columns with WT
  d1_2 <- df1[,grepl("_ALT",colnames(df1)) & !grepl("THRESH",colnames(df1))]
  
  # get the names of the columns where the SV is present
  matches2 <- colnames(d1_2)[which(d1_2>=1)]
  
  # count in how many tissue the SV is found in each mouse
  matches2_tx <- table(paste0("WT_",str_split_fixed(matches2,"_",3)[,2]))
  
  return(matches2_tx)
}
getsampleTx <- function(df1) {
  
  # get only columns with WT
  d1_2 <- df1[,grepl("_ALT",colnames(df1)) & !grepl("THRESH",colnames(df1))]
  
  # get the names of the columns where the SV is present
  matches2 <- colnames(d1_2)[which(d1_2>=1)]
  matches2 <- gsub("_ALT","",matches2)
  
  return(matches2)
}
annotateGenomicFeatures2 <- function(df,priorities,txdbGR,quiet,is.shuffle,txdb_raw) {
  
  # create SV id and convert file to GR
  df$id <- paste0(df$chrom,"-",df$start,"-",df$end,"-",df$sample_comparison)
  
  # get genomic features
  feature_Grl <- txdbGR
  dfgr2 <- annotateGR(df=df,TxDb=feature_Grl,quiet,prior=priorities,is.shuffle,txdb_raw)
  
  return(dfgr2)
}
wrapShuffleCount2 <- function(df,txdb,prior,txdb_GR,bootstrap,threads,sample_id,is.shuffle,txdb_raw) {
  counter <- 0
  
  for (sample in unique(df[,sample_id])) {
    
    counter <- counter + 1
    print(paste0("Analysing sample ",sample," (",counter,")"))
    
    df2 <- subset(df,df[,sample_id]==sample)
    
    # run function in parallel
    mergedl <- ( bplapply(seq_len(bootstrap),
                                   BPPARAM=MulticoreParam(worker=threads,RNGseed=seed_f),
                                   getExpectedGenomicFeatParallel2,
                                   df2,
                                   txdb,
                                   prior = c("Promoter","utr5","utr3","Exon","Intron","Downstream"),
                                   txdb_GR,
                                   is.shuffle=is.shuffle,
                                   txdb_raw=txdb_raw) )
    
    # merge
    for (i in 1:length(mergedl)) {
      if (i==1) {
        merged2 <- mergedl[[i]]
      } else {
        merged2 <- merge(merged2,mergedl[[i]],by="id",all=T)
      }
    }
    merged2[is.na(merged2)] <- 0
    
    # get only counts
    rownames(merged2) <- merged2$id
    merged2$id <- NULL

    # calc mean and sd
    nco <- ncol(merged2)
    merged2$Exp <- apply(merged2[,1:nco],1,mean)
    merged2$sd <- apply(merged2[,1:nco],1,sd)
    final <- merged2[,c("Exp","sd")]
    final$id <- rownames(merged2)
    
    # format as the obs file 
    final$sample_comparison <- str_split_fixed(final$id,"@",2)[,2]
    final$annotation <- str_split_fixed(final$id,"@",2)[,1]
    
    final <- final[,c("sample_comparison","annotation","Exp","sd")]
    final$bootstrap <- nco
    
    # bind
    if (counter==1) {
      f_final <- final
    } else {
      f_final <- rbind(f_final,final)
    }
  }
  
  return(f_final)
}
getExpectedGenomicFeatParallel2 <- function(i,df,txdb,prior,feature_Grl,is.shuffle,txdb_raw) {
  
  ## get shuffle
  #rand <- Shuf(bed = df,g = genome,bedtools = bedtools)
  #names(rand) <- names(df)
  ## bedtools is faster, but sometime, for unknown reasons, it generates files with some
  ## missing fields
  rand <- Shuf_ChIPseek(df,txdb)
  rand.ann <- annotateGenomicFeatures2(df=rand,
                                       priorities = prior,
                                       feature_Grl,
                                       quiet = T,
                                       is.shuffle,
                                       txdb_raw)
  exp <- as.data.frame(table(rand.ann[,c("annotation","sample_comparison")]))
  # get df with id and count only
  sample <- unique(exp$sample_comparison)
  exp$id <- paste0(exp$annotation,"@",exp$sample_comparison)
  exp <- exp[,c("id","Freq")]
  names(exp)[2] <- paste0("bootstrap",i)
  # check all the features are in output
  prior_df <- data.frame(id=paste0(prior,"@",sample))
  exp <- merge(prior_df,exp,by="id",all=T)
  exp <- replaceNAwithZero(exp)
  return(exp)
  
}
Shuf_ChIPseek <- function(bed,txdb) {
  # convert df to GR
  bedGR <- toGRanges(bed)
  # shuffle regions
  shufGR <- shuffle(bedGR,txdb)
  # convert GR back to df as the annotateGenomicFeatures expects df
  shuf <- as.data.frame(shufGR,row.names = NULL)
  # add missing columns (shuf only contains coord+strand)
  shuf <- shuf[,1:3]
  shuf <- cbind(shuf,bed[,4:ncol(bed)])
  # set colnames as input ones
  names(shuf) <- names(bed)
  return(shuf)
}
replaceNAwithZero <- function(dataf) {
  dataf[is.na(dataf)] <- 0
  return(dataf)
}
getOE <- function(real,shuffle) {
  
  # parse OBS
  obs <- as.data.frame(table(real[,c("annotation","sample_comparison")]))
  obs <- subset(obs,Freq>0)
  obs <- calcPerc2(obs,group = "sample_comparison",counts = "Freq")
  
  # merge obs exp
  merged <- merge(obs,shuffle,by='annotation')[,c("annotation","Freq","Exp","Tot","sd","bootstrap")]
  names(merged)[names(merged)=="Freq"] <- "Obs"
  
  # calc perc
  merged$Obs_perc <- merged$Obs / merged$Tot * 100
  merged$Exp_perc <- merged$Exp / merged$Tot * 100
  
  # calc zscore, pvalue, fdr
  merged <- zscoreCalc(merged)
  # format
  format1 <- c("Exp","sd","Obs_perc","Exp_perc","Enrichment","zscore")
  for (col in format1) {
    merged[,col] <- round(merged[,col],2)
  }
  format2 <- c("pvalue","fdr")
  for (col in format2) {
    merged[,col] <- formatC(merged[,col],digits = 2,format = 'e')
  }
  return(merged)
}
calcPerc2 <- function(df,group,counts) {
  df.aggr <- aggregate(df[,counts],by=list(df[,group]),sum)
  names(df.aggr) <- c(group,"Tot")
  df2 <- merge(df,df.aggr,by=group,all=T)
  if (!is.na(table(is.na(df2$Tot))["TRUE"])) { stop() }
  df2$Perc <- round(df2$Freq/df2$Tot*100,2)
  return(df2)
}
zscoreCalc <- function(dataf) {
  dataf$Enrichment <- log2 ( dataf$Obs / dataf$Exp )
  dataf$zscore <- (dataf$Obs - dataf$Exp) / dataf$sd
  dataf$zscore <- ifelse(is.na(dataf$zscore),0,dataf$zscore)
  dataf$pvalue <- 2*pnorm(-abs(dataf$zscore))
  dataf$fdr <- p.adjust(dataf$pvalue,method = 'fdr')
  return(dataf)
}
wrapShuffleCount_BND <- function(df,txdb,prior,txdb_GR,bootstrap,threads,sample_id,is.shuffle,txdb_raw) {
  
  counter <- 0
  
  for (sample in unique(df[,sample_id])) {
    
    counter <- counter + 1
    print(paste0("Analysing sample ",sample," (",counter,")"))
    
    df2 <- subset(df,df[,sample_id]==sample)
    
    # run function in parallel
    merged <- list.cbind( bplapply(seq_len(bootstrap),
                                   BPPARAM=MulticoreParam(worker=threads,RNGseed=seed_f),
                                   getExpectedGenomicFeatParallel_BND,
                                   df2,
                                   txdb,
                                   prior = c("Promoter","utr5","utr3","Exon","Intron","Downstream"),
                                   txdb_GR,
                                   is.shuffle=is.shuffle,
                                   txdb_raw=txdb_raw) )
    
    # do some checks
    col1 <- merged[,1]
    col2 <- merged[,3]
    if (!all(col1==col2)) { stop() }
    
    # get only counts
    odd_col <- seq_len(ncol(merged)) %% 2
    merged2 <- merged[ , odd_col == 0]
    
    # calc mean and sd
    nco <- ncol(merged2)
    merged2$Exp <- apply(merged2[,1:nco],1,mean)
    merged2$sd <- apply(merged2[,1:nco],1,sd)
    final <- merged2[,c("Exp","sd")]
    final$id <- col1
    
    # format as the obs file 
    final$sample_comparison <- str_split_fixed(final$id,"@",2)[,2]
    final$annotation <- str_split_fixed(final$id,"@",2)[,1]
    
    final <- final[,c("sample_comparison","annotation","Exp","sd")]
    final$bootstrap <- nco
    
    # bind
    if (counter==1) {
      f_final <- final
    } else {
      f_final <- rbind(f_final,final)
    }
  }
  
  return(f_final)
}
getExpectedGenomicFeatParallel_BND <- function(i,df,txdb,prior,feature_Grl,is.shuffle,txdb_raw) {
  
  ## get shuffle
  #rand <- Shuf(bed = df,g = genome,bedtools = bedtools)
  #names(rand) <- names(df)
  ## bedtools is faster, but sometime, for unknown reasons, it generates files with some
  ## missing fields
  rand <- Shuf_ChIPseek(df,txdb)
  rand.ann <- annotateGenomicFeaturesBND(df=rand,
                                         priorities = prior,
                                         feature_Grl,
                                         quiet = T,
                                         is.shuffle,
                                         txdb_raw)
  exp <- as.data.frame(table(rand.ann[,c("annotation","sample_comparison")]))
  # get df with id and count only
  sample <- unique(exp$sample_comparison)
  exp$id <- paste0(exp$annotation,"@",exp$sample_comparison)
  exp <- exp[,c("id","Freq")]
  names(exp)[2] <- paste0("bootstrap",i)
  # check all the features are in output
  prior2 <- c(prior,"Distal Intergenic")
  new_prior <- c()
  counter <- 0
  for (ii in 1:length(prior2)) {
    #print(ii)
    for (j in 1:length(prior2)) {
      counter <- counter + 1
      #print(j)
      if (prior2[[ii]] >= prior2[[j]]) {
        new_p <- paste0(prior2[[ii]],"/",prior2[[j]])
      } else {
        new_p <- paste0(prior2[[j]],"/",prior2[[ii]])
      }
      new_prior[counter] <- new_p
    }
  }
  new_prior <- unique(new_prior)
  prior_df <- data.frame(id=paste0(new_prior,"@",sample))
  exp <- merge(prior_df,exp,by="id",all=T)
  exp <- replaceNAwithZero(exp)
  return(exp)
  
}
wrap_countSRoverlap <- function(file,sr1,sr2,bedtools,i) {
  uno <- countSRoverlap(file,sr1,bedtools,i)
  names(uno)[2:3] <- paste0(names(uno)[2:3],"_1")
  # Simple rep dataset2
  due <- countSRoverlap(file,sr2,bedtools,i)
  names(due)[2:3] <- paste0(names(due)[2:3],"_2")
  system(paste0("rm ",file))
  return(list(uno,due))
}
wrap_GenomicFeat <- function(input,sv_common,sv_tissue) {
  
  ### COMMMON
  
  i <- 0
  for (dt in unique(input$somatic.anno)) {
    if (grepl("common",dt)) {
      for (svtype in sv_common) {
        i <- i + 1
        cat(paste0("\n\n0",i," - annotating ",dt," (",svtype,")\n"))
        
        # Get SVs of interest
        sv_of_i <- subset(input,somatic.anno==dt & SVTYPE==svtype)[,c("chrom","start","end","ID")]
        cat("\tselected",nrow(sv_of_i),svtype,"\n")
        
        ## annotate real SVs
        sv_of_i$sample_comparison <- dt
        sv_of_i.ann <- annotateGenomicFeatures2(df=sv_of_i,
                                                priorities = priorities,
                                                feature_Grl,
                                                quiet = T,
                                                is.shuffle = F,
                                                txdb_raw=txdb)
        # annotate shuffle SVs
        t1 <- Sys.time()
        exp <- wrapShuffleCount2(df = sv_of_i,
                                 txdb = txdb,
                                 prior = priorities,
                                 txdb_GR = feature_Grl,
                                 bootstrap = bootstrap,
                                 threads=5,
                                 sample_id="sample_comparison",
                                 is.shuffle = T,
                                 txdb_raw = txdb)
        t2 <- Sys.time()
        print(t2-t1)
        
        # combine in observed/expected table
        anno <- getOE(sv_of_i.ann,exp)
        anno$SVTYPE <- svtype
        anno$somatic.anno <- dt
        
        # bind
        if (i==1) {
          anno_common <- anno
        } else {
          anno_common <- rbind(anno_common,anno)
        }
        ## # write results
        ## write.table(anno,paste0(dt,"_OE_",svtype,"_",bootstrap,".txt"),row.names = F,sep='\t',quote=F)
      }
    }
  }
  
  
  #######################################################################################################################
  
  
  #######################################################################################################################
  
  ### TISSUE SPECIFIC (only BND)
  i <- i + 1
  dt <- "Somatic_tissue_private"
  sv_tissue <- "BND"
  
  cat(paste0("\n\n0",i," - annotating ",dt," (",sv_tissue,")\n"))
  
  # Get somatic common
  somatic_tissue <- subset(input,somatic.anno==dt)[,c("chrom","start","end","ID","SVTYPE","breakpoint2")]
  somatic_tissue <- subset(somatic_tissue,SVTYPE==sv_tissue)
  cat("\tselected",nrow(somatic_tissue),sv_tissue,"\n")
  
  # Divide by comparison
  comparison <- obj$tissue_specificity
  if (table(somatic_tissue$ID %in% comparison$SV_ID)["TRUE"]!=nrow(somatic_tissue)){stop()}
  somatic_tissue2 <- merge(comparison,somatic_tissue,by.x="SV_ID",by.y="ID",all.y=T)[,c("chrom","start","end","Type","breakpoint2")]
  print(table(somatic_tissue2$Type))
  
  
  ## annotate real SVs
  somatic_tissue2$sample_comparison <- somatic_tissue2$Type
  somatic_tissue2.ann <- annotateGenomicFeaturesBND(df=somatic_tissue2,
                                                    priorities = priorities,
                                                    feature_Grl,
                                                    quiet = T,
                                                    is.shuffle = F,
                                                    txdb_raw=txdb)
  # annotate shuffle SVs
  t1 <- Sys.time()
  exp <- wrapShuffleCount_BND(df = somatic_tissue2,
                              txdb = txdb,
                              prior = priorities,
                              txdb_GR = feature_Grl,
                              bootstrap = bootstrap,
                              threads=detectCores(),
                              sample_id="sample_comparison",
                              is.shuffle = T,
                              txdb_raw = txdb)
  t2 <- Sys.time()
  print(t2-t1)
  
  
  
  # Modify annotation by pasting genomic annotation annd tissue-specificity sample
  # this is done in order to safely run the function below getOE
  somatic_tissue2.ann$annotation <- paste0(somatic_tissue2.ann$sample_comparison,"@",somatic_tissue2.ann$annotation)
  exp$annotation <- paste0(exp$sample_comparison,"@",exp$annotation)
  
  # Get OE table
  somatic_tissue_OE <- getOE(somatic_tissue2.ann,exp)
  somatic_tissue_OE$SVTYPE <- sv_tissue
  somatic_tissue_OE$somatic.anno <- dt
  
  # Split anno and tissue-spec
  somatic_tissue_OE$type <- str_split_fixed(somatic_tissue_OE$annotation,"@",2)[,1]
  somatic_tissue_OE$annotation <- str_split_fixed(somatic_tissue_OE$annotation,"@",2)[,2]
  
  
  ## write.table(somatic_tissue_OE,paste0(dt,"_OE_BND_",bootstrap,".txt"),row.names = F,sep='\t',quote=F)
  
  # bind common and tissue-specific
  anno_common$type <- "common"
  ff <- rbind(anno_common,somatic_tissue_OE)
  return(ff)
}
wrapGO_tissue_specific <- function(background,svtype,svclass,fdr,gp_version) {
  
  
  ### Background - all the genes that host a genic somatic variant
  cat("\nGetting background genes (all the genes that host a genic somatic variant)")
  
  # All somatic variants
  somatic <- background
  
  # Gene name relative to genes hosting somatic variants
  all_genes_genic_variant1 <- unique(subset(somatic,annotation1!="Distal Intergenic" &
                                              annotation1!="Downstream")[,"gene_name1"])
  all_genes_genic_variant2 <- unique(subset(somatic,annotation2!="Distal Intergenic" &
                                              annotation2!="Downstream")[,"gene_name2"])
  all_genes_genic_variant <- unique(c(all_genes_genic_variant1,all_genes_genic_variant2))
  # discard empty (not all the variants have a 2nd gene, but only the BND)
  all_genes_genic_variant <- all_genes_genic_variant[!all_genes_genic_variant==""]
  cat(": ",length(all_genes_genic_variant),"\n")
  
  
  
  
  
  
  ### TISSUE SPECIFIC (only BND)
  i <- 1
  
  dt <- svclass
  
  cat(paste0("\n\n0",i," - GO on ",dt," (",svtype,")\n"))
  
  # Get somatic tissue specific private, only BND
  somatic_tissue <- subset(somatic,somatic.anno==dt)
  somatic_tissue <- subset(somatic_tissue,SVTYPE==svtype)
  cat("\tselected",nrow(somatic_tissue),svtype,"\n")
  
  # Divide by pairwise comparisons tissue by tissue
  comparison <- obj$tissue_specificity
  if (table(somatic_tissue$ID %in% comparison$SV_ID)["TRUE"]!=nrow(somatic_tissue)){stop()}
  somatic_tissue2 <- merge(comparison,somatic_tissue,by.x="SV_ID",by.y="ID",all.y=T)
  #print(table(somatic_tissue2$Type))
  
  
  ###### Start GO
  
  # Get columns of interest (annotation and gene name)
  input_df <- unique(somatic_tissue2[,c("Type","annotation1","annotation2","gene_name1","gene_name2")])
  
  
  # Comparison by comparison GO analysis
  features <- c("Promoter","Exon","Intron")
  counter <- 0
  for (tp in unique(input_df$Type)) {
    cat("\n\t",tp,"\n")
    
    # select variants of the comparison of interest
    d1 <- subset(input_df,Type==tp)
    
    for (feat in features) {
      cat("\t\t",feat,":")
      
      # Get genes associated to the feature
      d1_genes_1 <- unique(subset(d1,annotation1==feat)[,"gene_name1"])
      d1_genes_2 <- unique(subset(d1,annotation2==feat)[,"gene_name2"])
      d1_genes <- unique(c(d1_genes_1,d1_genes_2))
      cat(length(d1_genes),"\n")
      # run go
      n_significant_genes <- round(5*length(d1_genes)/100,0)
      cat("\t\tGO must be found in at least 5% of the genes (",n_significant_genes,"genes)\n")
      cat("\t\tFDR =",fdr,"\n")
      res <- doGprof2(query_gene = d1_genes,bg_genes = all_genes_genic_variant,
                      fdr = fdr,GO_name = "GO:BP",n_significant_genes = n_significant_genes,gp_version = gp_version)
      # bind
      if (nrow(res)>0) {
        counter <- counter + 1
        res$Type <- tp
        res$annotation <- feat
        if (counter==1) {
          final_res <- res
        } else {
          final_res <- rbind(final_res,res)
        }
      }
    }
  }
  return(final_res)
}
wrapGO_common <- function(background,svtype,svclass,fdr,gp_version) {
  
  ### Background - all the genes that host a genic somatic variant
  cat("\nGetting background genes (all the genes that host a genic somatic variant)")
  
  # All somatic variants
  somatic <- background
  
  # Gene name relative to genes hosting somatic variants
  all_genes_genic_variant1 <- unique(subset(somatic,annotation1!="Distal Intergenic" &
                                              annotation1!="Downstream")[,"gene_name1"])
  all_genes_genic_variant2 <- unique(subset(somatic,annotation2!="Distal Intergenic" &
                                              annotation2!="Downstream")[,"gene_name2"])
  all_genes_genic_variant <- unique(c(all_genes_genic_variant1,all_genes_genic_variant2))
  # discard empty (not all the variants have a 2nd gene, but only the BND)
  all_genes_genic_variant <- all_genes_genic_variant[!all_genes_genic_variant==""]
  cat(": ",length(all_genes_genic_variant),"\n")
  
  
  ### START ANALYSIS
  
  dt <- svclass
  
  cat(paste0("\n\n- GO on ",dt," (",svtype,")\n"))
  
  # Get variants of interest
  somatic_tissue <- subset(somatic,somatic.anno==dt)
  somatic_tissue <- subset(somatic_tissue,SVTYPE==svtype)
  cat("\tselected",nrow(somatic_tissue),svtype,"\n")
  
  ###### Start GO
  
  # Get columns of interest (annotation and gene name)
  input_df <- unique(somatic_tissue[,c("ID","annotation1","gene_name1")])
  
  
  # Comparison by comparison GO analysis
  features <- c("Promoter","Exon","Intron")
  counter <- 0
  for (feat in features) {
    cat("\t\t",feat,":")
    
    # Get genes associated to the feature
    d1_genes <- unique(subset(input_df,annotation1==feat)[,"gene_name1"])
    cat(length(d1_genes),"\n")
    # run go
    n_significant_genes <- round(5*length(d1_genes)/100,0)
    cat("\t\tGO must be found in at least 5% of the genes (",n_significant_genes,"genes)\n")
    cat("\t\tFDR =",fdr,"\n")
    res <- doGprof2(query_gene = d1_genes,bg_genes = all_genes_genic_variant,
                    fdr = fdr,GO_name = "GO:BP",n_significant_genes = n_significant_genes,gp_version = gp_version)
    # bind
    if (nrow(res)>0) {
      counter <- counter + 1
      res$annotation <- feat
      if (counter==1) {
        final_res <- res
      } else {
        final_res <- rbind(final_res,res)
      }
    }
  }
  if (!exists("final_res")) {final_res=data.frame()}
  return(final_res)
}
annotateGenomicFeaturesBND2 <- function(df,priorities,txdbGR,quiet,is.shuffle,txdb_raw) {
  
  # create SV id and convert file to GR
  df$id <- paste0(df$chrom,"-",df$start,"-",df$end,"-",df$sample_comparison)
  
  # get genomic features
  feature_Grl <- txdbGR
  dfgr2 <- annotateGR_BND2(df=df,TxDb=feature_Grl,quiet,prior=priorities,is.shuffle,txdb_raw)
  
  return(dfgr2)
}
annotateGR_BND2 <- function(df,TxDb,quiet,prior,is.shuffle,txdb_raw) {
  ### Define priorities
  if (length(prior) != 6) {
    stop("the vector of priorities does not contain 6 features")
  }
  if (! quiet) {
    print(paste0("Assigning  priorities as follow: ",prior[1],", ",prior[2],", ",prior[3],", ",prior[4],", ",
                 prior[5],", ",prior[6]))
  }
  
  ######## Breakpoint 1
  
  ### Convert to granges
  gr_bed_1 <- toGRanges(df)
  ### Overlap gr_bed with feature_Grl - return is dataframe (chr,start,end,id,txname,annotation)
  overlaps_1 <- grOverlapTxDb(GR_bed=gr_bed_1,grl_txdb=TxDb,c_prior=prior,quiet=quiet,is.shuffle=is.shuffle,txdb=txdb_raw)
  
  
  ######## Breakpoint 2
  
  ### Convert to granges
  df2 <- df
  df2$chrom <- str_split_fixed(df2$breakpoint2,":",2)[,1]
  df2$start <- as.numeric(str_split_fixed(df2$breakpoint2,":",2)[,2])
  df2$end <- df2$start + 1
  
  gr_bed_2 <- toGRanges(df2)
  ### Overlap gr_bed with feature_Grl - return is dataframe (chr,start,end,id,txname,annotation)
  overlaps_2 <- grOverlapTxDb(GR_bed=gr_bed_2,grl_txdb=TxDb,c_prior=prior,quiet=quiet,is.shuffle=is.shuffle,txdb=txdb_raw)
  
  
  ######## Merge Breakpoint 1 and 2
  if ( table(overlaps_1$id %in% overlaps_2$id)["TRUE"]!=nrow(overlaps_1)){ stop() }
  overlaps <- merge(overlaps_1,overlaps_2,by="id")
  overlaps$annotation <- paste0(overlaps$annotation.x,"/",overlaps$annotation.y)
  
  
  ### Add the annotation to the original file
  if (is.shuffle) {
    overlaps <- overlaps[,c("id","annotation")]
  } else {
    #overlaps_1 <- overlaps_1[,c("id","tx_name2","tx_coord","distanceTSS","distanceTx","annotation")]
    #overlaps_2 <- overlaps_2[,c("id","tx_name2","tx_coord","distanceTSS","distanceTx","annotation")]
    # unique annotation formatting
    overlaps$tx_name2 <- paste0(overlaps$tx_name2.x,"/",overlaps$tx_name2.y)
    overlaps$tx_coord <- paste0(overlaps$tx_coord.x,"/",overlaps$tx_coord.y)
    overlaps$distanceTSS <- paste0(overlaps$distanceTSS.x,"/",overlaps$distanceTSS.y)
    overlaps$distanceTx <- paste0(overlaps$distanceTx.x,"/",overlaps$distanceTx.y)
    overlaps <- overlaps[,c("id","annotation","tx_name2","tx_coord","distanceTSS","distanceTx")]
    colnames(overlaps)[colnames(overlaps)=="tx_name2"] <- "tx_name"
  }
  final <- merge(df,overlaps,by="id",sort=F)
  
  
  ### Check input and output have the same length
  if (nrow(final)!=nrow(df)) {
    stop("input and output do not have the same length")
  }
  return(final)
}
getGeneSymbol <- function(df,anno) {
  # check all the tx in vcf are also in annotation file
  if (!table(df$tx_name %in% anno$tx_name)["TRUE"]==nrow(df)) { 
    stop("some transcripts in vcf are missing from annotation file") 
  }
  out <- merge(df,anno,by="tx_name",all.x=T,sort=F)
  return(out)
}
mergeBND <- function(orig,bnd1,bnd2,orig_sort) {
  # merge original df containing BND with df containing gene name associated to breakp 1 and breakp 2
  mm <- merge(orig,merge(bnd1[,c("ID","gene_name")],
                         bnd2[,c("ID","gene_name")],all=T,by="ID",sort=F),
              all=T,by="ID",sort=F)
  mm$gene_name <- paste0(mm$gene_name.x,"/",mm$gene_name.y)
  mm$gene_name.x=NULL
  mm$gene_name.y=NULL
  # order col as original
  mm <- mm[,c(colnames(orig_sort),"gene_name")]
  return(mm)
}
mergeBND_SR <- function(sr1_res,sr2_res,bp) {
  # merge BND breakpoint 1 and breakpoint 2 - TRF SR
  rr <- merge(sr1_res,sr2_res,by='ID',all=T,sort=F)
  # SR name
  rr[,paste0("SR_",bp,".x")] <- ifelse(is.na(rr[,paste0("SR_",bp,".x")]),"No SR",rr[,paste0("SR_",bp,".x")])
  rr[,paste0("SR_",bp,".y")] <- ifelse(is.na(rr[,paste0("SR_",bp,".y")]),"No SR",rr[,paste0("SR_",bp,".y")])
  rr[,paste0("SR_",bp)] <- paste0(rr[,paste0("SR_",bp,".x")],"/",rr[,paste0("SR_",bp,".y")])
  # SR coords
  rr[,paste0("SR_ID",".x")] <- ifelse(is.na(rr[,paste0("SR_ID",".x")]),"No SR",rr[,paste0("SR_ID",".x")])
  rr[,paste0("SR_ID",".y")] <- ifelse(is.na(rr[,paste0("SR_ID",".y")]),"No SR",rr[,paste0("SR_ID",".y")])
  rr[,paste0("SR_ID",bp)] <- paste0(rr[,paste0("SR_ID",".x")],"/",rr[,paste0("SR_ID",".y")])
  # intersect is 1, this is just 1 bp
  rr[,paste0("intersect_",bp)] <- 1
  rr <- rr[,c("ID",paste0("SR_",bp),paste0("intersect_",bp),paste0("SR_ID",bp))]
  return(rr)
}
getCountsBND_SR <- function(dataf,group,item_to_count) {
  # count
  ct <- as.data.frame(table(dataf[,c(group,item_to_count)]))
  # calc %
  ct <- calcPerc2(ct,group,"Freq")
  return(ct)
}
getShuffleBND_SR_counts2 <- function(n,dataf,path_genome,path_sr_shuf,path_simple_rep1,path_simple_rep2,bedtools,bootstrap) {
  
  ### ***START***
  
  # shuffle coords
  seed_i <- 6492+n
  system(paste0("bedtools shuffle -g ",path_genome," -seed ",seed_i," -incl ",path_sr_shuf," -i uno.tmp > uno.shuf_boot",n,".tmp"))
  
  # intersect with SR
  test <- wrap_countSRoverlap(file = paste0("uno.shuf_boot",n,".tmp"),sr1 = path_simple_rep1,sr2 = path_simple_rep2,bedtools,i=n)
  test1 <- test[[1]] # overlap with TRF SR
  test2 <- test[[2]] # overlap with rmsk SR
  # merge TRF SR and rmsk SR
  test_int <- merge(test1,test2,by="ID",all=T)
  if (nrow(test_int)!=nrow(rep)) {stop()}
  # Define whether selecting first or second SR annotation - priority to TRF
  test_int$SR <- ifelse(is.na(test_int$SR_1),test_int$SR_2, # if first is NA, get second
                        ifelse(is.na(test_int$SR_2),test_int$SR_1, # if second is NA, get first
                               ifelse(test_int$SR_1!=test_int$SR_2,test_int$SR_1, # if first is different from second, get first
                                      ifelse(test_int$SR_1==test_int$SR_2,test_int$SR_1,"ERROR"))) )
  # add SR motif len
  test_int$SR_motif_len <- nchar(test_int$SR) - 3 # -3 as we dont want to count "(", ")" and "n"
  
  # get cols of interest ID2-->Type and SR_motif_len
  test_int$Type <- str_split_fixed(test_int$ID,"@",3)[,1]
  test_int <- test_int[,c("Type","SR_motif_len")]
  # count
  exp <- getCountsBND_SR(dataf=test_int,group="Type",item_to_count="SR_motif_len")
  exp$SR_motif_len <- paste0("n=",exp$SR_motif_len)
  exp <- exp[,c("SR_motif_len","Freq","Type",'Tot')]
  exp$id <- paste0(exp$SR_motif_len,"@",exp$Type)
  exp <- exp[,c("id","Freq","Tot")]
  names(exp)[2] <- paste0("boot_",n)
  names(exp)[3] <- paste0("Tot_boot_",n)
  
  return(exp)
}
bindShuf <- function(inp_ls,obs_df,name) {
  
  for (i in 1:length(inp_ls)) {
    
    if (i%%100==0) { print(paste0("Bootstrap:",i))}
    
    # getdata
    df <- inp_ls[[i]]
    # merge
    if (i!=1) {
      ff <- merge(ff,df,by="id",all=T,sort=F)
      # check totals are equal
      if (all(ff$Tot_boot_1!=ff[,paste0("Tot_boot_",i)],na.rm = T)) {stop()}
      # remove column with totals - just leave one
      ff[,paste0("Tot_boot_",i)] <- NULL
    } else {
      ff <- df
    }
  }
  
  # calc mean and sd
  ff$Exp <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,mean)
  ff$sd <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,sd)
  ff2 <- ff[,c("id","Tot_boot_1","Exp","sd")]
  
  
  # bind with observed
  obs_df$id <- paste0("n=",obs_df[,name],"@",obs_df$Type)
  merged <- merge(obs_df,ff2,all.x=T,sort=F,by="id")
  # check toals obs and exp are equal
  if (all(merged$Tot==merged$Tot_boot_1,na.rm = T)) { merged$Tot_boot_1=NULL } else { stop() }
  
  # ge motifs found at least once
  merged <- subset(merged,Freq>0)
  
  # convert NA to 0
  merged <- replaceNAwithZero(merged)
  colnames(merged)[colnames(merged)=="Freq"] <- "Obs"
  
  # calc zscore
  merged <- zscoreCalc(merged)
  
  return(merged)
  
}
bindShuf_motif <- function(inp_ls,obs_df,name) {
  
  for (i in 1:length(inp_ls)) {
    
    if (i%%100==0) { print(paste0("Bootstrap:",i))}
    
    # getdata
    df <- inp_ls[[i]]
    # merge
    if (i!=1) {
      ff <- merge(ff,df,by="id",all=T,sort=F)
      # check totals are equal
      if (all(ff$Tot_boot_1!=ff[,paste0("Tot_boot_",i)],na.rm = T)) {stop()}
      # remove column with totals - just leave one
      ff[,paste0("Tot_boot_",i)] <- NULL
    } else {
      ff <- df
    }
  }
  
  # calc mean and sd
  ff$Exp <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,mean)
  ff$sd <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,sd)
  ff2 <- ff[,c("id","Tot_boot_1","Exp","sd")]
  
  
  # bind with observed
  obs_df$id <- paste0(obs_df[,name],"@",obs_df$Type)
  merged <- merge(obs_df,ff2,all.x=T,sort=F,by="id")
  # check toals obs and exp are equal
  if (all(merged$Tot==merged$Tot_boot_1,na.rm = T)) { merged$Tot_boot_1=NULL } else { stop() }
  
  # ge motifs found at least once
  merged <- subset(merged,Freq>0)
  
  # convert NA to 0
  merged <- replaceNAwithZero(merged)
  colnames(merged)[colnames(merged)=="Freq"] <- "Obs"
  
  # calc zscore
  merged <- zscoreCalc(merged)
  
  return(merged)
  
}
aggreCounts <- function(df) {
  ### aggregate counts
  # Obs %
  df2 <- aggregate(df$Perc,by=list(df$Type,df$SR_len_class),sum)
  names(df2) <- c("Type","SR_len_class","Obs_perc")
  # Exp %
  df2$Exp_perc <- aggregate(df$Exp_perc,by=list(df$Type,df$SR_len_class),sum)[,"x"]
  # Obs
  df2$Obs <- aggregate(df$Obs,by=list(df$Type,df$SR_len_class),sum)[,"x"]
  # Exp
  df2$Exp <- aggregate(df$Exp,by=list(df$Type,df$SR_len_class),sum)[,"x"]
  # sd / var %
  df$var_perc <- df$sd_perc ^ 2
  df2$var_perc <- aggregate(df$var_perc,by=list(df$Type,df$SR_len_class),sum)[,"x"]
  df2$sd_perc <- sqrt(df2$var_perc)
  # sd / var
  df$var <- df$sd ^ 2
  df2$var <- aggregate(df$var,by=list(df$Type,df$SR_len_class),sum)[,"x"]
  df2$sd <- sqrt(df2$var)
  # enrichment log2
  df2$Enrichment <- log2 ( df2$Obs_perc / df2$Exp_perc )
  # zscore
  df2$zscore <- (df2$Obs_perc - df2$Exp_perc) / df2$sd_perc
  df2$zscore <- ifelse(is.na(df2$zscore),0,df2$zscore)
  # pval
  df2$pvalue <- 2*pnorm(-abs(df2$zscore))
  # fdr
  df2$fdr <- p.adjust(df2$pvalue,method = 'fdr')
  
  return(df2)
}
getShuffleBND_SR_counts2_SR2 <- function(n,dataf,path_genome,path_sr_shuf,path_simple_rep1,path_simple_rep2,bedtools,bootstrap) {
  
  ### ***START***
  
  ###################### *** SHUFFFLING ***
  # shuffle coords
  seed_i <- 6492+n
  system(paste0("bedtools shuffle -g ",path_genome," -seed ",seed_i," -incl sr.bed -i uno.tmp > uno.shuf_boot",n,".tmp"))
  ###################### *** END ***
  
  
  ###################### *** INTERSECTING WITH SR ***
  
  # intersect with SR
  test <- wrap_countSRoverlap(file = paste0("uno.shuf_boot",n,".tmp"),
                              sr1 = path_simple_rep1,
                              sr2 = path_simple_rep2,
                              bedtools,i=n)
  test1 <- test[[1]] # overlap with TRF SR
  test2 <- test[[2]] # overlap with rmsk SR
  # merge TRF SR and rmsk SR
  test_int <- merge(test1,test2,by="ID",all=T)
  if (nrow(test_int)!=nrow(dataf)) {stop()}
  # Define whether selecting first or second SR annotation - priority to TRF
  test_int$SR <- ifelse(is.na(test_int$SR_1),test_int$SR_2, # if first is NA, get second
                        ifelse(is.na(test_int$SR_2),test_int$SR_1, # if second is NA, get first
                               ifelse(test_int$SR_1!=test_int$SR_2,test_int$SR_1, # if first is different from second, get first
                                      ifelse(test_int$SR_1==test_int$SR_2,test_int$SR_1,"ERROR"))) )
  
  # get cols of interest ID2-->Type and SR_motif_len
  test_int$Type <- str_split_fixed(test_int$ID,"@",3)[,1]
  test_int <- test_int[,c("Type","SR")]
  # count
  exp <- getCountsBND_SR(dataf=test_int,group="Type",item_to_count="SR")
  exp$id <- paste0(exp$SR,"@",exp$Type)
  exp <- exp[,c("id","Freq","Tot")]
  names(exp)[2] <- paste0("boot_",n)
  names(exp)[3] <- paste0("Tot_boot_",n)
  return(exp)
}
getShuffleCommonDEL_SR_counts2_SR2 <- function(n,dataf,path_genome,path_sr_shuf,path_simple_rep1,path_simple_rep2,bedtools,bootstrap) {
  
  ### ***START***
  
  ###################### *** SHUFFFLING ***
  # shuffle coords
  seed_i <- 6492+n
  system(paste0("bedtools shuffle -g ",path_genome," -seed ",seed_i," -incl sr.bed -i uno.tmp > uno.shuf_boot",n,".tmp"))
  ###################### *** END ***
  
  
  ###################### *** INTERSECTING WITH SR ***
  
  # intersect with SR
  test <- wrap_countSRoverlap_shuf_common(file = paste0("uno.shuf_boot",n,".tmp"),
                                          sr1 = path_simple_rep1,
                                          sr2 = path_simple_rep2,
                                          bedtools,i=n)
  test1 <- test[[1]] # overlap with TRF SR
  test2 <- test[[2]] # overlap with rmsk SR
  # merge TRF SR and rmsk SR
  test_int <- merge(test1,test2,by="ID",all=T)
  if (nrow(test_int)!=nrow(dataf)) {stop()}
  # Define whether selecting first or second SR annotation - priority to TRF
  test_int$SR <- ifelse(is.na(test_int$SR_1),test_int$SR_2, # if first is NA, get second
                        ifelse(is.na(test_int$SR_2),test_int$SR_1, # if second is NA, get first
                               ifelse(test_int$SR_1!=test_int$SR_2,test_int$SR_1, # if first is different from second, get first
                                      ifelse(test_int$SR_1==test_int$SR_2,test_int$SR_1,"ERROR"))) )
  
  # get cols of interest ID2-->Type and SR
  test_int$somatic.anno <- str_split_fixed(test_int$ID,"@",2)[,1]
  test_int <- test_int[,c("somatic.anno","SR")]
  # count
  exp <- getCountsBND_SR(dataf=test_int,group="somatic.anno",item_to_count="SR")
  exp$id <- paste0(exp$SR,"@",exp$somatic.anno)
  exp <- exp[,c("id","Freq","Tot")]
  names(exp)[2] <- paste0("boot_",n)
  names(exp)[3] <- paste0("Tot_boot_",n)
  return(exp)
}
getNN_SR <- function(path_to_SR,idx) {
  # subset SR to get only (NN)n
  sr1bed <- fread(path_to_SR,data.table = F)
  sr1bed$len <- nchar(sr1bed$V4) -3
  sr1bed <- subset(sr1bed,len==2)
  sr1bed$len <- NULL
  fwrite(sr1bed,paste0("sr",idx,".bed.tmp"),na = NA,row.names = F,col.names = F,sep='\t',quote=F)
  return(paste0("sr",idx,".bed.tmp"))
}
countSRoverlap_shuf_common <- function(path_dataf,rmsk,bedtools,i) {
  
  # intersect SV with TE coords
  cmd <- paste0(bedtools," intersect -wo -a ",path_dataf," -b ",rmsk," > tmp.intersect_bootstr",i)
  #cat("Executing",cmd,"\n")
  system(cmd)
  # read results
  intersect.orig <- fread(paste0("tmp.intersect_bootstr",i),header=F,sep='\t',data.table = F)
  system(paste0("rm tmp.intersect_bootstr",i))
  
  # assign column names
  #col1 <- colnames(fread(path_dataf,header=T,sep='\t'))
  col1 <- c("chrom","start","end","ID")
  col2 <- c("chr_SR","start_SR","end_SR","SR","SR2","strand_SR","intersect")
  header1 <- c(col1,col2)
  if (length(header1)!=ncol(intersect.orig)) {stop("error in the column formatting")}
  colnames(intersect.orig) <- header1
  
  # get len and ID
  intersect.orig$len <- intersect.orig$end_SR-intersect.orig$start_SR
  intersect.orig$SR_ID <- paste(intersect.orig$chr_SR,intersect.orig$start_SR,intersect.orig$end_SR,sep='-')
  
  # select columns of interest
  intersect <- unique(intersect.orig[,c("ID","SR","intersect","len","SR_ID")])
  # remove redundancies (one SV can overlap more than one TE)
  selected <- as.data.frame(intersect %>% 
                              group_by(ID) %>% 
                              slice_max(order_by = tibble(intersect,len), n = 1,with_ties=F) ) # first select the one with the longest intersect, if intersection length is the same, get the variant overlapping the longest SR
  if (nrow(selected)!=length(unique(intersect$ID))) { stop("Redundant SVs found") }
  # Add info to melted file
  selected <- selected[,c("ID","SR","intersect","SR_ID")]
  
  return(selected)
}
wrap_countSRoverlap_shuf_common <- function(file,sr1,sr2,bedtools,i) {
  uno <- countSRoverlap_shuf_common(file,sr1,bedtools,i)
  names(uno)[2:3] <- paste0(names(uno)[2:3],"_1")
  # Simple rep dataset2
  due <- countSRoverlap_shuf_common(file,sr2,bedtools,i)
  names(due)[2:3] <- paste0(names(due)[2:3],"_2")
  system(paste0("rm ",file))
  return(list(uno,due))
}
getShuffleBND_SR_counts2_common <- function(n,dataf,path_genome,path_sr_shuf,path_simple_rep1,path_simple_rep2,bedtools,bootstrap) {
  
  
  ### ***START***
  
  # shuffle coords
  seed_i <- 6492+n
  system(paste0("bedtools shuffle -g ",path_genome," -seed ",seed_i," -incl ",path_sr_shuf," -i uno.tmp > uno.shuf_boot",n,".tmp"))
  
  # intersect with SR
  test <- wrap_countSRoverlap_shuf_common(file = paste0("uno.shuf_boot",n,".tmp"),sr1 = path_simple_rep1,sr2 = path_simple_rep2,bedtools,i=n)
  test1 <- test[[1]] # overlap with TRF SR
  test2 <- test[[2]] # overlap with rmsk SR
  # merge TRF SR and rmsk SR
  test_int <- merge(test1,test2,by="ID",all=T)
  if (nrow(test_int)!=nrow(common)) {stop()}
  # Define whether selecting first or second SR annotation - priority to TRF
  test_int$SR <- ifelse(is.na(test_int$SR_1),test_int$SR_2, # if first is NA, get second
                        ifelse(is.na(test_int$SR_2),test_int$SR_1, # if second is NA, get first
                               ifelse(test_int$SR_1!=test_int$SR_2,test_int$SR_1, # if first is different from second, get first
                                      ifelse(test_int$SR_1==test_int$SR_2,test_int$SR_1,"ERROR"))) )
  # add SR motif len
  test_int$SR_motif_len <- nchar(test_int$SR) - 3 # -3 as we dont want to count "(", ")" and "n"
  
  # get cols of interest ID2-->Type and SR_motif_len
  test_int$somatic.anno <- str_split_fixed(test_int$ID,"@",3)[,1]
  test_int <- test_int[,c("somatic.anno","SR_motif_len")]
  # count
  exp <- getCountsBND_SR(dataf=test_int,group="somatic.anno",item_to_count="SR_motif_len")
  exp$SR_motif_len <- paste0("n=",exp$SR_motif_len)
  exp <- exp[,c("SR_motif_len","Freq","somatic.anno",'Tot')]
  exp$id <- paste0(exp$SR_motif_len,"@",exp$somatic.anno)
  exp <- exp[,c("id","Freq","Tot")]
  names(exp)[2] <- paste0("boot_",n)
  names(exp)[3] <- paste0("Tot_boot_",n)
  
  return(exp)
}
bindShuf_common <- function(inp_ls,obs_df,name) {
  
  for (i in 1:length(inp_ls)) {
    
    if (i%%100==0) { print(paste0("Bootstrap:",i))}
    
    # getdata
    df <- inp_ls[[i]]
    # merge
    if (i!=1) {
      ff <- merge(ff,df,by="id",all=T,sort=F)
      # check totals are equal
      if (all(ff$Tot_boot_1!=ff[,paste0("Tot_boot_",i)],na.rm = T)) {stop()}
      # remove column with totals - just leave one
      ff[,paste0("Tot_boot_",i)] <- NULL
    } else {
      ff <- df
    }
  }
  
  # calc mean and sd
  ff$Exp <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,mean)
  ff$sd <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,sd)
  ff2 <- ff[,c("id","Tot_boot_1","Exp","sd")]
  
  
  # bind with observed
  obs_df$id <- paste0("n=",obs_df[,name],"@",obs_df$somatic.anno)
  merged <- merge(obs_df,ff2,all.x=T,sort=F,by="id")
  # check toals obs and exp are equal
  if (all(merged$Tot==merged$Tot_boot_1,na.rm = T)) { merged$Tot_boot_1=NULL } else { stop() }
  
  # ge motifs found at least once
  merged <- subset(merged,Freq>0)
  
  # convert NA to 0
  merged <- replaceNAwithZero(merged)
  colnames(merged)[colnames(merged)=="Freq"] <- "Obs"
  
  # calc zscore
  merged <- zscoreCalc(merged)
  
  return(merged)
  
}
aggreCounts2 <- function(df) {
  ### aggregate counts
  # Obs %
  df2 <- aggregate(df$Perc,by=list(df$somatic.anno,df$SR_len_class),sum)
  names(df2) <- c("somatic.anno","SR_len_class","Obs_perc")
  # Exp %
  df2$Exp_perc <- aggregate(df$Exp_perc,by=list(df$somatic.anno,df$SR_len_class),sum)[,"x"]
  # Obs
  df2$Obs <- aggregate(df$Obs,by=list(df$somatic.anno,df$SR_len_class),sum)[,"x"]
  # Exp
  df2$Exp <- aggregate(df$Exp,by=list(df$somatic.anno,df$SR_len_class),sum)[,"x"]
  # sd / var %
  df$var_perc <- df$sd_perc ^ 2
  df2$var_perc <- aggregate(df$var_perc,by=list(df$somatic.anno,df$SR_len_class),sum)[,"x"]
  df2$sd_perc <- sqrt(df2$var_perc)
  # sd / var
  df$var <- df$sd ^ 2
  df2$var <- aggregate(df$var,by=list(df$somatic.anno,df$SR_len_class),sum)[,"x"]
  df2$sd <- sqrt(df2$var)
  # enrichment log2
  df2$Enrichment <- log2 ( df2$Obs_perc / df2$Exp_perc )
  # zscore
  df2$zscore <- (df2$Obs_perc - df2$Exp_perc) / df2$sd_perc
  df2$zscore <- ifelse(is.na(df2$zscore),0,df2$zscore)
  # pval
  df2$pvalue <- 2*pnorm(-abs(df2$zscore))
  # fdr
  df2$fdr <- p.adjust(df2$pvalue,method = 'fdr')
  
  return(df2)
}
bindShuf_motif_common <- function(inp_ls,obs_df,name) {
  
  for (i in 1:length(inp_ls)) {
    
    if (i%%100==0) { print(paste0("Bootstrap:",i))}
    
    # getdata
    df <- inp_ls[[i]]
    # merge
    if (i!=1) {
      ff <- merge(ff,df,by="id",all=T,sort=F)
      # check totals are equal
      if (all(ff$Tot_boot_1!=ff[,paste0("Tot_boot_",i)],na.rm = T)) {stop()}
      # remove column with totals - just leave one
      ff[,paste0("Tot_boot_",i)] <- NULL
    } else {
      ff <- df
    }
  }
  
  # calc mean and sd
  ff$Exp <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,mean)
  ff$sd <- apply(ff[,!grepl("Tot_boot",colnames(ff)) & grepl("boot_",colnames(ff))],1,sd)
  ff2 <- ff[,c("id","Tot_boot_1","Exp","sd")]
  
  
  # bind with observed
  obs_df$id <- paste0(obs_df[,name],"@",obs_df$somatic.anno)
  merged <- merge(obs_df,ff2,all.x=T,sort=F,by="id")
  # check toals obs and exp are equal
  if (all(merged$Tot==merged$Tot_boot_1,na.rm = T)) { merged$Tot_boot_1=NULL } else { stop() }
  
  # ge motifs found at least once
  merged <- subset(merged,Freq>0)
  
  # convert NA to 0
  merged <- replaceNAwithZero(merged)
  colnames(merged)[colnames(merged)=="Freq"] <- "Obs"
  
  # calc zscore
  merged <- zscoreCalc(merged)
  
  return(merged)
  
}
getN_SR <- function(path_to_SR,idx) {
  # subset SR to get only (NN)n
  sr1bed <- fread(path_to_SR,data.table = F)
  sr1bed$len <- nchar(sr1bed$V4) -3
  sr1bed <- subset(sr1bed,len==1)
  sr1bed$len <- NULL
  fwrite(sr1bed,paste0("sr",idx,".bed.tmp"),na = NA,row.names = F,col.names = F,sep='\t',quote=F)
  return(paste0("sr",idx,".bed.tmp"))
}
wrap_overlapBND_TEs <- function(bnd_bed,path_to_rmsk,bedtools,i,TEsplit) {
  
  # overlap BND and TEs
  teo <- overlapBND_TEs(bnd_bed,path_to_rmsk,bedtools,i)
  
  # count
  count_teo <- count_overlapBND_TEs(teo,TEsplit)
  
  return(count_teo)
}
overlapBND_TEs <- function(bnd_bed,path_to_rmsk,bedtools,i) {
  
  # write bed to file
  write.table(bnd_bed,paste0("cicci_",i),row.names = F,col.names = F,sep='\t',quote=F)
  path_dataf <- paste0("cicci_",i)
  
  # intersect SV with TE coords
  cmd <- paste0(bedtools," intersect -wao -a ",path_dataf," -b ",path_to_rmsk," > tmp.intersect_bootstr",i)
  system(cmd)
  # read results
  intersect.orig <- fread(paste0("tmp.intersect_bootstr",i),header=F,sep='\t',data.table = F)
  system(paste0("rm cicci_",i," tmp.intersect_bootstr",i))
  
  # assign column names
  col1 <- colnames(bnd_bed)
  col2 <- c("chr_TE","start_TE","end_TE","TE1","TE2","strand_TE","intersect")
  header1 <- c(col1,col2)
  if (length(header1)!=ncol(intersect.orig)) {stop("error in the column formatting")}
  colnames(intersect.orig) <- header1
  
  # select columns of interest
  intersect.orig$TElen <- intersect.orig$end_TE - intersect.orig$start_TE
  intersect <- unique(intersect.orig[,c("ID","type","Type","TE1","TE2","intersect","TElen")])
  intersect$ID <- paste0(intersect$ID,"_",intersect$type,"_",intersect$Type)
  intersect$type <- NULL
  intersect$Type <- NULL
  intersect$TE <- paste0(intersect$TE1,":",intersect$TE2)
  names(intersect) <- c("ID","discard1","discard2","overlap","TElen","TE")
  intersect$discard1 <- NULL
  intersect$discard2 <- NULL
  # assign -1 to SV not overlapping TE
  intersect$TE <- ifelse(intersect$TE==".:-1","No TE",intersect$TE)
  # remove redundancies (one BND can overlap more than one TE)
  # since the BND are all long 1bp, in case of redundancies, select the BND overlapping the longest TEs
  selected <- as.data.frame(intersect %>% 
                              group_by(ID) %>% 
                              slice_max(TElen, n = 1,with_ties=F) )
  if (nrow(selected)!=length(unique(intersect$ID))) { stop("Redundant SVs found") }
  # Add info to melted file
  selected <- selected[,c("ID","TE")]
  
  return(selected)
  
}
count_overlapBND_TEs <- function(overlap_df,TEsplit) {
  
  ### split TE annotation
  overlap_df$TE <- ifelse(overlap_df$TE=="No TE","No TE:No TE:No TE",overlap_df$TE)
  overlap_df$TE <- str_split_fixed(overlap_df$TE,":",3)[,TEsplit]
  
  ### merge breakp 1 and 2
  # bp1
  bp1 <- overlap_df[grepl("_bp1",overlap_df$ID),]
  bp1$ID <- gsub("bp1_","",bp1$ID)
  # bp1
  bp2 <- overlap_df[grepl("_bp2",overlap_df$ID),]
  bp2$ID <- gsub("bp2_","",bp2$ID)
  # merge
  merged <- merge(bp1,bp2,by="ID")
  if (nrow(merged)!=nrow(bp1)) {stop()}
  if (nrow(merged)!=nrow(bp2)) {stop()}
  # sorted annotation
  merged$TE <- ifelse(merged$TE.x>merged$TE.y,paste0(merged$TE.x,"/",merged$TE.y),paste0(merged$TE.y,"/",merged$TE.x))
  
  # count by sample
  merged$sample <- str_split_fixed(merged$ID,"_",2)[,2]
  count_bnd <- as.data.frame(table(merged[,c("TE","sample")]))
  
  return(count_bnd)
  
}
wrap_overlapBND_TEs_shuf <- function(bnd_bed,txdb,path_to_rmsk,bedtools,TEsplit,bootstrap) {
  # run function in parallel
  merged_res <- bplapply(seq_len(bootstrap),
                         BPPARAM=MulticoreParam(worker=multicoreWorkers(),RNGseed=seed_f),
                         getExpectedoverlapBND_TEs,
                         bnd_bed,
                         txdb,path_to_rmsk,bedtools,
                         TEsplit)
  
  # merge all the bootstrap res
  for (i in 1:length(merged_res)) {
    if (i==1) {
      exp_res <- merged_res[[i]]
    } else {
      exp_res <- merge(exp_res,merged_res[[i]],by="ID",all=T,sort=F)
    }
  }
  # calc mean and sd
  nco <- ncol(exp_res)
  exp_res$Expected <- apply(exp_res[,2:nco],1,mean)
  exp_res$sd <- apply(exp_res[,2:nco],1,sd)
  exp_res2 <- exp_res[,c("ID","Expected","sd")]
  exp_res2$bootstrap <- nco-1
  
  exp_res2$TE <- str_split_fixed(exp_res2$ID,"@",2)[,1]
  exp_res2$sample <- str_split_fixed(exp_res2$ID,"@",2)[,2]
  exp_res2$ID <- NULL
  
  return(exp_res2)
}
getExpectedoverlapBND_TEs <- function(i,df,txdb,path_to_rmsk,bedtools,TEsplit) {
  
  if (i%%10==0) {
    cat(paste0("\tBootstrap: ",i,"\n"))
  }
  
  ## get shuffle
  #rand <- Shuf(bed = df,g = genome,bedtools = bedtools)
  #names(rand) <- names(df)
  ## bedtools is faster, but sometime, for unknown reasons, it generates files with some
  ## missing fields
  rand <- Shuf_ChIPseek(df,txdb)
  
  # overlap TEs + count overlap
  exp <- wrap_overlapBND_TEs(rand,path_to_rmsk,bedtools,i,TEsplit=2)
  exp$ID <- paste0(exp$TE,"@",exp$sample)
  names(exp)[3] <- paste0("bootstrap",i)
  
  exp$TE <- NULL
  exp$sample <- NULL
  
  return(exp)
  
}
addComparisonType <- function(dataf,df_comparison) {
  
  # Add comparison type
  dataf2 <- merge(dataf,df_comparison,by.x="ID",by.y="SV_ID",all.x=T,sort=F)
  
  # Expand dataframe based on the "Samples" column (which contains a comma-separated list of the samples where the sv is found)
  dataf2_expand <- as.data.frame(tidyr::separate_rows(dataf2, Samples, sep=","))
  
  # Discard samples not matching comparison type
  # there might be cases where the samples are KID_6,CTX_2 and the comparisons KID-vs-HIP and CTX-vs-HIP
  # in such cases, when the df is expanded, we need to select the matching samples and comparisons
  dataf2_expand$tissue_sample <- str_split_fixed(dataf2_expand$Samples,"_",2)[,1]
  dataf2_expand$tissue_type <- str_split_fixed(dataf2_expand$Type,"-",3)[,1]
  dataf2_expand <- subset(dataf2_expand,tissue_sample==tissue_type)
  
  # Define non-redundant comparison type
  dataf2_expand$Type2 <- 
    ifelse(dataf2_expand$Type=="HIP-vs-KID" |dataf2_expand$Type=="KID-vs-HIP","HIP vs\nKID",
           ifelse(dataf2_expand$Type=="CTX-vs-KID" |dataf2_expand$Type=="KID-vs-CTX","CTX vs\nKID",
                  ifelse(dataf2_expand$Type=="HIP-vs-CTX" |dataf2_expand$Type=="CTX-vs-HIP","CTX vs\nHIP","error")))
  
  # Modify names in a nicer way ...
  dataf2_expand$somatic.anno <- ifelse(dataf2_expand$somatic.anno=="Somatic_tissue_private","Somatic private",
                                       ifelse(dataf2_expand$somatic.anno=="Somatic_tissue_recurrent","Somatic recurrent","error"))
  return(dataf2_expand)
}
convertToBed <- function(dataf) {
  # breakpoint 1
  dataf_1 <- dataf
  dataf_1$chrom <- str_split_fixed(str_split_fixed(dataf_1$SV_ID,"@",2)[,1],":",2)[,1]
  dataf_1$start <- as.numeric(str_split_fixed(str_split_fixed(dataf_1$SV_ID,"@",2)[,1],":",2)[,2])
  dataf_1$end   <- dataf_1$start + 1
  
  # breakpoint 2
  dataf_2 <- dataf
  dataf_2$chrom <- str_split_fixed(str_split_fixed(dataf_2$SV_ID,"@",2)[,2],":",2)[,1]
  dataf_2$start <- as.numeric(str_split_fixed(str_split_fixed(dataf_2$SV_ID,"@",2)[,2],":",3)[,2])
  dataf_2$end   <- dataf_2$start + 1
  
  # bed file
  bnd_bed <- rbind(dataf_1,dataf_2)[,c("chrom","start","end","SV_ID")]
  
  return(bnd_bed)
}
bedtools_closest <- function(df1,df2) {
  # write to file
  write.table(df1,"bedA.bed",row.names = F,col.names = F,sep='\t',quote=F)
  write.table(df2,"bedB.bed",row.names = F,col.names = F,sep='\t',quote=F)
  # sort
  system("bedtools sort -i bedA.bed > bedA.bedS")
  system("bedtools sort -i bedB.bed > bedB.bedS")
  # closest
  system("bedtools closest -d -t all -nonamecheck -a bedA.bedS -b bedB.bedS > closest.bed")
  clos <- read.table("closest.bed",header = F,sep='\t')
  names(clos) <- c(paste0(colnames(df1),"_ont"),paste0(colnames(df2),"_illumina"),"distance")
  # clean
  system("rm bedA.bed bedB.bed bedA.bedS bedB.bedS closest.bed")
  return(clos)
  
}
sortChromID <- function(svid) {
  chr1 <- str_split_fixed(str_split_fixed(svid,"@",2)[,1],":",2)[,1]
  chr2 <- str_split_fixed(str_split_fixed(svid,"@",2)[,2],":",2)[,1]
  chr  <- ifelse(chr1>chr2,paste0(chr1,"@",chr2),paste0(chr2,"@",chr1))
  return(chr)
}
validatingBND <- function(ont,ill,D) {
  
  # format
  colnames(ont)[grepl("SV_ID",colnames(ont))] <- "SV_ID"
  colnames(ill)[grepl("SV_ID",colnames(ill))] <- "SV_ID"
  
  # get bnd
  bnd_ont <- subset(ont,SVTYPE_ont=="BND")
  ill_ont <- ill[grepl("BND",ill$SVTYPE),]
  
  # convert to bed
  bnd_ont_bed <- convertToBed(bnd_ont)
  ill_ont_bed <- convertToBed(ill_ont)
  
  # get closest
  closest <- unique(bedtools_closest(bnd_ont_bed,ill_ont_bed))
  
  # threshold on D
  closest <- subset(closest,distance<D)
  
  # in case of dupl, get the bnd with the same chromosome combination
  closest$chr_mOnt <- sortChromID(closest$SV_ID_ont)
  closest$chr_mIll <- sortChromID(closest$SV_ID_illumina)
  closest$equal.chr <- ifelse(closest$chr_mOnt==closest$chr_mIll,1,0)
  
  closest <- as.data.frame(closest %>% 
                             group_by(SV_ID_ont) %>% 
                             slice_max(equal.chr, n = 1,with_ties=T) )
  
  # if still dupl, get the closest one
  closest <- as.data.frame(closest %>% 
                             group_by(SV_ID_ont) %>% 
                             slice_min(distance, n = 1,with_ties=F) )
  
  # get ID of the validated ones
  validated <- unique(closest[,c("SV_ID_ont","SV_ID_illumina")])
  
  # print msg
  cat(paste0("\tvalidation %: ",round(nrow(validated)/nrow(bnd_ont)*100),"% (",nrow(validated),"/",nrow(bnd_ont),")\n"))
  
  return(validated)
}
getQueryGenes <- function(background,svtype,svclass) {
  
  
  ### Background - all the genes that host a genic somatic variant
  
  # All somatic variants
  somatic <- background
  
  # Gene name relative to genes hosting somatic variants
  all_genes_genic_variant1 <- unique(subset(somatic,annotation1!="Distal Intergenic" &
                                              annotation1!="Downstream")[,"gene_name1"])
  all_genes_genic_variant2 <- unique(subset(somatic,annotation2!="Distal Intergenic" &
                                              annotation2!="Downstream")[,"gene_name2"])
  all_genes_genic_variant <- unique(c(all_genes_genic_variant1,all_genes_genic_variant2))
  # discard empty (not all the variants have a 2nd gene, but only the BND)
  all_genes_genic_variant <- all_genes_genic_variant[!all_genes_genic_variant==""]
  
  ### TISSUE SPECIFIC (only BND)
  i <- 1
  
  dt <- svclass
  
  # Get somatic tissue specific private, only BND
  somatic_tissue <- subset(somatic,somatic.anno==dt)
  somatic_tissue <- subset(somatic_tissue,SVTYPE==svtype)
  
  # Divide by pairwise comparisons tissue by tissue
  comparison <- obj$tissue_specificity
  if (table(somatic_tissue$ID %in% comparison$SV_ID)["TRUE"]!=nrow(somatic_tissue)){stop()}
  somatic_tissue2 <- merge(comparison,somatic_tissue,by.x="SV_ID",by.y="ID",all.y=T)
  
  
  ###### Start GO
  
  # Get columns of interest (annotation and gene name)
  input_df <- unique(somatic_tissue2[,c("Type","annotation1","annotation2","gene_name1","gene_name2")])
  
  
  # Comparison by comparison GO analysis
  features <- c("Promoter","Exon","Intron")
  counter <- 0
  res_go  <- data.frame()
  for (tp in unique(input_df$Type)) {

    # select variants of the comparison of interest
    d1 <- subset(input_df,Type==tp)
    
    for (feat in features) {
      
      # Get genes associated to the feature
      d1_genes_1 <- unique(subset(d1,annotation1==feat)[,"gene_name1"])
      d1_genes_2 <- unique(subset(d1,annotation2==feat)[,"gene_name2"])
      d1_genes <- unique(c(d1_genes_1,d1_genes_2))
      
      # save to df
      df_go <- data.frame(feat_go = feat,
                          Type    = tp,
                          genes   = d1_genes)
      res_go <- rbind(res_go,df_go)
      
    }
  }
  
  # add background genes
  bg_genes <- data.frame(feat_go = "background",
                         Type    = "background",
                         genes   = all_genes_genic_variant)
  res_go2 <- rbind(bg_genes,res_go)
  return(res_go2)
}
getQueryGenes_common <- function(background,svtype,svclass) {
  
  
  ### Background - all the genes that host a genic somatic variant
  
  # All somatic variants
  somatic <- background
  
  # Gene name relative to genes hosting somatic variants
  all_genes_genic_variant1 <- unique(subset(somatic,annotation1!="Distal Intergenic" &
                                              annotation1!="Downstream")[,"gene_name1"])
  all_genes_genic_variant2 <- unique(subset(somatic,annotation2!="Distal Intergenic" &
                                              annotation2!="Downstream")[,"gene_name2"])
  all_genes_genic_variant <- unique(c(all_genes_genic_variant1,all_genes_genic_variant2))
  # discard empty (not all the variants have a 2nd gene, but only the BND)
  all_genes_genic_variant <- all_genes_genic_variant[!all_genes_genic_variant==""]
  
  ### TISSUE SPECIFIC (only BND)
  i <- 1
  
  dt <- svclass
  
  # Get somatic Somatic_common_private DEL
  somatic_tissue <- subset(somatic,somatic.anno==dt)
  somatic_tissue <- subset(somatic_tissue,SVTYPE==svtype)
  

  
  # Get columns of interest (annotation and gene name)
  input_df <- unique(somatic_tissue[,c("ID","annotation1","gene_name1")])
  
  
  # Comparison by comparison GO analysis
  features <- c("Promoter","Exon","Intron")
  counter <- 0
  res_go  <- data.frame()
  for (feat in features) {
    
    # Get genes associated to the feature
    d1_genes <- unique(subset(input_df,annotation1==feat)[,"gene_name1"])
    
    # save to df
    df_go <- data.frame(feat_go      = feat,
                        SVTYPE       = svtype,
                        somatic.anno = dt,
                        genes        = d1_genes)
    res_go <- rbind(res_go,df_go)
    
    
  }
  
  # add background genes
  bg_genes <- data.frame(feat_go      = "background",
                         SVTYPE       = svtype,
                         somatic.anno = dt,
                         genes        = all_genes_genic_variant)
  res_go2 <- rbind(bg_genes,res_go)
  return(res_go2)
}
getGOanno <- function(bgenes,genome,geneID) {
  
  
  # get GO annotation
  genome <- "mm39"        # only selects org.Mm.eg.db for the GO lookup
  geneID     <- "geneSymbol"
  # fetched once rather than re-fetched inside each of the 9 goseq() calls
  gene2cat <- getgo(bgenes, genome, geneID, fetch.cats = "GO:BP")
  names(gene2cat) <- tolower(names(gene2cat))
  
  # getgo() marks non-BP terms "Other" and only drops them when a gene has more
  # than one entry, so a gene annotated solely in CC/MF keeps a lone "Other".
  # GO:0008150 is the BP root, reached by ancestor propagation, and informative
  # for nothing.
  gene2cat <- lapply(gene2cat, function(x) {
    if (is.null(x)) return(NULL)
    x <- x[!is.na(x) & x != "Other" & x != "GO:0008150"]
    if (!length(x)) NULL else x
  })
  gene2cat <- gene2cat[!vapply(gene2cat, is.null, logical(1))]
  
  # the reverse --> GO associated to genes
  cat2gene <- split(rep(names(gene2cat), lengths(gene2cat)),
                    unlist(gene2cat, use.names = FALSE))
  
  out_ls <- list(gene2cat,cat2gene)
  names(out_ls) <- c("gene_to_cat","cat_to_gene")
  return(out_ls)
}
getGeneLength <- function(gtf_genes) {
  gene_l_anno             <- unique(gtf_genes[,c("gene_name","width")])
  gene_l_anno$gene_name   <- tolower(gene_l_anno$gene_name)
  gene_l_anno_ag          <- aggregate(gene_l_anno$width,by=list(gene_l_anno$gene_name),max)
  names(gene_l_anno_ag)   <- c("gene_name","width")
  return(gene_l_anno_ag)
}
runGOseq <- function(genes_oi,bggenes,gene_l,genecat,catgene,min_genes) {

  # 1 = query gene, 0 = rest of the background
  # all the d1 genes are in the universe, so we can directly classify as 0/1 the universe genes
  de <- as.integer(bggenes %in% genes_oi$genes)
  names(de) <- bggenes
  if (sum(de) != length(unique(genes_oi$genes))) stop("query genes missing from background")
  
  # get gene length
  bias_len <- subset(gene_l,gene_name %in% names(de))
  de <- de[bias_len$gene_name]
  if (all(bias_len$gene_name==names(de))) {
    bias_len_c <- bias_len$width
    names(bias_len_c) <- names(de)
  } else {
    stop()
  }
  
  # Fitting the Probability Weighting Function (PWF)
  pwf <- nullp(de, bias.data = bias_len_c, plot.fit = F)
  
  # corrected GO analysis
  corrected <- goseq(pwf, gene2cat = genecat, method = "Wallenius")
  #corrected <- goseq(pwf, genome = "mm39",id = "ensGene", method = "Wallenius")
  #corrected2 <- goseq(pwf, genome = "mm39",id = "geneSymbol", method = "Wallenius")
  
  # remove genes supported by <5 genes
  corrected_f <- subset(corrected,numDEInCat>=min_genes)
  
  # get info of interest
  rr     <- corrected_f[, c("category","over_represented_pvalue","numDEInCat","numInCat","term")]
  rr$fdr <- p.adjust(rr$over_represented_pvalue, "BH")
  rr     <- rr[order(rr$fdr), ]
  
  # calc enrichment
  N <- nrow(pwf)        # genes in the universe
  K <- sum(pwf$DEgenes) # query genes among them
  rr$Enrichment <- (rr$numDEInCat / K) / (rr$numInCat / N)
  
  
  # get genes associated to the GO
  de_genes  <- names(de)[de == 1]
  rr$genes <- vapply(rr$category,
                     function(k) paste(sort(unname(sym_map[intersect(catgene[[k]], de_genes)])),
                                       collapse = ","),
                     character(1))
  return(rr)
}
