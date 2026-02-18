library(fsusieR)
plot(out$res)
out$res$cs
plot(out$res$lBF[[1]])
  mean(out$start_bin- out$end_bin)

  
  out$Y
  
  fsusie_log_plot(obj=out$res,
                  chr = out$chr,
                  log1p_count = TRUE,
                  pos0 = out$pos[1],
                  pos1 = max(out$pos),
                  X=out$X,
                  Y=out$Y,
                  snp_info=out$info_SNP, cs =1  ) 
  
  
  
  
  
  
  library(AnnotationHub)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(Gviz)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  ah <- AnnotationHub()
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  # This plot visualizes the effect of a chosen credible set (CS) using
  # the Gviz Bioconductor package. The input "obj" is the fsusie object
  # (i.e., the result of running fsusie).  X and Y are the genotype and
  # (scaled) read count matrices, respectively.
  fsusie_log_plot <- function (obj, chr, pos0, pos1, X, Y, snp_info, cs = 1,
                               log1p_count=FALSE,
                               
                               data_splice=NULL,
                               plot_cred_band=TRUE,
                               type_data="p",positions) {
    
    # Extract the relevant genes and exons in the specified region
    region_genes <- genes(txdb,columns = c("tx_id","gene_id"))
    
    # Subset the genes and exons to the region of interest.
    region_genes <- subsetByOverlaps(region_genes,
                                     GRanges(seqnames = chr,
                                             ranges = IRanges(pos0,pos1)))
    
    # Generate a sequence of positions with a length of 1,024.
    if (missing (positions)){
      
      positions <- seq(pos0,pos1,length.out = 1024)
      
    }
    markers <- obj$cs[[cs]]
    j       <- which.max(obj$pip[markers])
    marker  <- markers[j]
    x       <- X[,marker]
    
    
    
    
    if(log1p_count){
      if (! length(which(x==2))>1){
        read_counts <- rbind(colMeans(log1p(Y[x == 0,])),
                             colMeans(log1p(Y[x == 1,])) )
      }else{
        read_counts <- rbind(colMeans(log1p(Y[x == 0,])),
                             colMeans(log1p(Y[x == 1,])),
                             colMeans(log1p(Y[x == 2,])))
      }
      
    }else{
      if (! length(which(x==2))>1){
        
        read_counts <- rbind(colMeans(Y[x == 0,]),
                             colMeans(Y[x == 1,]) )
      }else{
        read_counts <- rbind(colMeans(Y[x == 0,]),
                             colMeans(Y[x == 1,]),
                             colMeans(Y[x == 2,]))
      }
      
    }
    
    
    if (plot_cred_band){
      
      effect=rbind (obj$fitted_func[[cs]],
                    obj$cred_band[[cs ]],
                    rep(0, length(obj$fitted_func[[cs]])))
      
      group_cred= c(1:3,0)
      group_colors <- c("black" ,"royalblue","royalblue","royalblue" )
    }else{
      effect=obj$fitted_func[[cs]]
      
      
    }
    
    
    
    
    
    # Create a "data track" to show the CS effect.
    cex <- 0.6
    
    
    
    effect_track <-
      DataTrack(range = GRanges(seqnames = chr,
                                ranges = IRanges(start = positions,
                                                 end = positions + 1)),
                data = effect, genome = "hg38",
                groups= group_cred,
                name = paste("CS",cs),type = "l",col = group_colors,
                track.margin = 0.05,cex.title = cex,cex.axis = cex,
                col.axis = "black",col.title = "black",
                fontface = "plain",background.title = "white",
                fontface.title = 1)
    
    
    
    
    # Create another "data track" to show the read counts.
    
    n0  <- sum(x == 0)
    n1  <- sum(x == 1)
    n2  <- sum(x == 2)
    id  <- snp_info[marker,"ID"]
    ref <- snp_info[marker,"REF"]
    alt <- snp_info[marker,"ALT"]
    
    if (! length(which(x==2))>1){
      groups <- c(sprintf("%s %s%s (n = %d)",id,ref,ref,n0),
                  sprintf("%s %s%s (n = %d)",id,ref,alt,n1) )
      geno_colors <- c("navyblue","turquoise" )
    }else{
      groups <- c(sprintf("%s %s%s (n = %d)",id,ref,ref,n0),
                  sprintf("%s %s%s (n = %d)",id,ref,alt,n1),
                  sprintf("%s %s%s (n = %d)",id,alt,alt,n2))
      geno_colors <- c("navyblue","turquoise","darkorange")
    }
    
    
    if (mean(effect) > 0) {
      groups <- factor(groups,rev(groups))
      geno_colors <- rev(geno_colors)
    } else {
      groups <- factor(groups,groups)
    }
    
    type_data="p"
    lab_y =ifelse(log1p_count, "avg. log1p count","avg. count")
    temp_pos = seq(pos0,pos1,length.out =dim(read_counts)[2])
    data_track <- DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = temp_pos ,
                                                             end =  temp_pos + 1)),
                            data = read_counts,genome = "hg38",
                            groups = groups,
                            
                            name = lab_y  , type = type_data, #"p",#type = "l",
                            col = geno_colors  ,
                            track.margin = 0.05,cex.title = cex,cex.axis = cex,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",background.title = "white",
                            fontface.title = 1,cex.legend = cex, cex=0.2)
    
    # Create an "ideogram" track.
    ideo_track <- IdeogramTrack(genome = "hg38",chromosome = chr)
    
    # Create a "genome axis" track.
    genome_track <- GenomeAxisTrack(col.axis = "black",col.title = "black")
    
    # Create a "gene region" track.
    gene_track <- GeneRegionTrack(txdb,genome = "hg38",chromosome = chr,
                                  pos0 = pos0,pos1 = pos1,name = "",
                                  showId = TRUE,geneSymbol = TRUE,
                                  col.axis = "black",col.title = "black",
                                  transcriptAnnotation = "symbol",
                                  rotation.title = 0,cex.title = cex,
                                  col = "salmon",fill = "salmon",
                                  background.title = "white")
    
    # Map gene IDs to gene symbols.
    gene_ids <- unique(unlist(region_genes$gene_id))
    
    # Map to gene symbols using org.Hs.eg.db
    gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,keys = gene_ids,
                                          columns = "SYMBOL",
                                          keytype = "ENTREZID")
    n <- nrow(gene_symbols)
    if (n > 0) {
      for (i in 1:n) {
        j <- which(gene_track@range@elementMetadata@listData$gene ==
                     gene_symbols$ENTREZID[i])
        gene_track@range@elementMetadata@listData$id[j]     <- gene_symbols$SYMBOL[i]
        gene_track@range@elementMetadata@listData$symbol[j] <- gene_symbols$SYMBOL[i]
      }
    }
    
    
    if( !is.null(data_splice)){
      junction_ranges <- GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = as.numeric(sub(".*_(\\d+)_.*", "\\1", data_splice$Name)),
          end = as.numeric(sub(".*_(\\d+)$", "\\1", data_splice$Name))
        ),
        names = data_splice$Description
      )
      
      # Offset overlapping junctions for better visualization
      # Create AnnotationTrack for splicing junctions with adjusted thickness
      junction_track <- AnnotationTrack(
        range = junction_ranges,
        genome = "hg38",
        chromosome = chr,
        name = "Splicing Junctions",
        stacking = "squish",
        col = "darkgreen",
        fill = "lightgreen",
        track.margin = 0.05,
        cex.title = cex,
        cex.axis = cex,
        col.axis = "black",
        col.title = "black",
        fontface = "plain",
        background.title = "white",
        height = 0.3  # Adjust height here to make rectangles thinner
      )
      # Combine all tracks into a single plot
      tracks <- c(
        ideo_track,
        genome_track,
        effect_track,
        data_track,
        gene_track,
        junction_track
      )
      
      # Plot the tracks
      return(plotTracks(tracks, from = pos0, to = pos1, sizes = c(1, 1.75, 2, 4, 5, 1)))
    }else{
      # Combine all tracks into a single plot.
      tracks <- c(ideo_track,
                  genome_track,
                  effect_track,
                  data_track,
                  gene_track)
      return(plotTracks(tracks,from = pos0,to = pos1,sizes = c(1,1.75,2,4,5)))
    }
    
    
  }
  