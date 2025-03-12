path="D:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie"
source(paste0(path,"/code/plot_all_effect_log.R"))
source(paste0(path,"/code/plot_log.R"))

load(paste0(path,"/output/local/ENSG00000080824.csv.gz.RData"))
plot_fsusie_log (out = out, log_effect = FALSE   ) 


library(fsusieR)

plot_susiF_pip(out$res, pos_SNP = as.numeric(out$info_SNP$POS))+
  
  geom_rect(aes(xmin = out$locus[1], 
                xmax =out$locus[2],
                ymin = -0.01,
                ymax = 1.02), 
            alpha = 0.0, color = "red")


#load(paste0(path,"/output/df_brain_cortex_junction.RData")) 
#fsusie_log_plot(out$res,chr = paste0("chr",out$chr),
#                pos0 = out$locus[1],pos1 = out$locus[2],
#                out$X,out$Y,snp_info = out$info_SNP,cs = 1,
#                effect_log=TRUE,
#                log1p_count=TRUE )

#fsusie_log_plot(out$res,chr = paste0("chr",out$chr),
#                pos0 = out$locus[1],pos1 = out$locus[2],
#                out$X,out$Y,snp_info = out$info_SNP,cs = 2,
#                effect_log=TRUE,
#                log1p_count=TRUE )
 







obj= out$res
Y=log1p(as.matrix(out$Y/out$size_factor_local) )
X=as.matrix(out$X)
to="TI"
verbose=TRUE
max_scale=10
filter_cs =FALSE
filter.number = 10
family =  "DaubLeAsymm"


post_processing=to
if (is.null( obj$pos))
{
  pos <- 1:dim(Y)[2]
}else{
  pos=obj$pos
}

names_colX <-  colnames(X)
tidx <- which(apply(X,2,var)==0)
if( length(tidx)>0){
  warning(paste("Some of the columns of X are constants, we removed" ,length(tidx), "columns"))
  X <- X[,-tidx]
}
map_data <- remap_data(Y=Y,
                       pos= pos,
                       verbose=vebose,
                       max_scale=max_scale)

outing_grid <- map_data$outing_grid
Y           <- map_data$Y
X <- colScale(X)

indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))
# centering input
#Y0 <-  colScale(Y , scale=FALSE)
Y  <- colScale(Y )
# out <- out_prep(     obj            = obj,
#
#                     X             = X,
#                     indx_lst      = indx_lst,
#                     filter_cs     = filter_cs,
#                      outing_grid   = outing_grid,
#                     filter.number = filter.number,
#                     family        = family,
#                      post_processing=  post_processing,
#                     tidx          = tidx,
#                     names_colX    = names_colX,
#                     pos           = pos
#)

obj1 <-   smash_regression(
  obj,
  Y             =     sweep(Y  , 2, attr(Y , "scaled:scale"),  "*"),
  X             = X,
  alpah=0.01
)



plot( obj1$fitted_func[[1]])

out1=out

obj=obj1
chr= paste0("chr",out1$chr)
pos0 =out1$locus[1]
pos1=out1$locus[2]
X=out1$X 
Y=out1$Y
snp_info=out1$info_SNP
cs = 1
log1p_count=TRUE
data_splice=NULL
plot_cred_band=TRUE
type_data="p"

# Extract the relevant genes and exons in the specified region
region_genes <- genes(txdb,columns = c("tx_id","gene_id"))

# Subset the genes and exons to the region of interest.
region_genes <- subsetByOverlaps(region_genes,
                                 GRanges(seqnames = chr,
                                         ranges = IRanges(pos0,pos1)))

# Generate a sequence of positions with a length of 1,024.
positions <- seq(pos0,pos1,length.out = 1024)

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

### CS1 -----

uni_res= univariate_functional_regression(Y=log1p(Y ),
                                          X=as.matrix(X[,obj$cs[[cs]][1]],
                                                      ncol=1),
                                          method="TI"
)
if (plot_cred_band){
  
  effect=rbind (uni_res$effect_estimate ,
                uni_res$cred_band ,
                rep(0, length(obj$fitted_func[[cs]])))
  
  group_cred= c(1:3,0)
  group_colors <- c("black" ,"royalblue","royalblue","royalblue" )
}else{
  effect=obj$fitted_func[[cs]]
  
  
}



# Create a "data track" to show the CS effect.
cex <- 1

group_lwd= c(1,2,1,1)
group_lty= c(1,1,2,2)
effect_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data = effect, genome = "hg38",
            groups= group_cred,
            lwd=group_lwd,
            lty=group_lty,
            name = paste("Effect CS",cs),type = "l",col = group_colors,
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )


plotTracks(effect_track)
# Create another "data track" to show the read counts.



n0  <- sum(x == 0)
n1  <- sum(x == 1)
n2  <- sum(x == 2)
id  <- "rs376197939" #snp_info[marker,"ID"]
ref <- "AA"#snp_info[marker,"REF"]
alt <- "AT" #snp_info[marker,"ALT"]

if (! length(which(x==2))>1){
  groups <- c(sprintf("\t    %s %s (n = %d)",id, ref,n0),
              sprintf("\t    %s %s (n = %d)",id, alt,n1) )
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


lab_y =ifelse(log1p_count, "Avg. log1p count","Avg. count")
data_track <- DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions,
                                                         end = positions + 1)),
                        data = read_counts,genome = "hg38",
                        groups = groups,
                        name = lab_y  , type = type_data, #"p",#type = "l",
                        col = geno_colors  ,
                        track.margin = 0.05,cex.title = cex,cex.axis = cex,
                        
                        col.axis = "black",col.title = "black",
                        fontface = "plain",background.title = "white",
                        fontface.title = 1,cex.legend = cex, cex=2)

data_track <-DataTrack(range = GRanges(seqnames = chr,
                                       ranges = IRanges(start = positions,
                                                        end = positions + 1)),
                       data = read_counts,genome = "hg38",
                       groups = groups,
                       name = lab_y , type = type_data,  col = geno_colors ,track.margin = 0.05,cex.title = cex,cex.axis = cex,
                       
                       col.axis = "black",col.title = "black",
                       fontface = "plain",background.title = "white", 
                       fontface.title = 1, cex= .6,cex.legend = 1.1)


plotTracks(data_track)

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
                              rotation.title = 0,cex.title = 2,
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


# Combine all tracks into a single plot.
tracks <- c(ideo_track,
            genome_track,
            effect_track,
            data_track,
            gene_track)


plotTracks(tracks,from = pos0,to = pos1,sizes = c(1,1.75,2,4,2)) 

folder_path=  paste0(getwd(),
                     "/plot/"
)
file_path <- file.path(folder_path, "HSP90AA1_cs1.pdf")
pdf(file_path, width =11.69, height = 8.27 )  # A4 in inches



plotTracks(tracks,from = pos0-200,
           to = pos1+200,
           sizes = c(1,1 ,4,4,2),
           cex.main=1.2, cex.title = 1.
) 
dev.off()




