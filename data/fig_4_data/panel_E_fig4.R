path= getwd()
data = readRDS(paste0(path , 
                      "/data/fig_4_data/Fig4_data.rds"))


# Interface: Assign datasets dynamically
plot_df <- data$e[[1]]
haQTL_df <- data$e[[2]]
gene_info <- data$e[[3]]
sumstat <- data$e[[4]]
pip_df <- data$e[[5]]

# Custom Labeller
custom_labeller <- function(x) {
  x %>% 
    gsub("DeJager_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", " \n ", .)             
}

# Parameters
text_size <- 20
view_win <- c(207317782, 207895513)
gene <- c("chr1:205117782-208795513", "1_205972031_208461272", "ENSG00000203710","ENSG00000117322")

# -----------------
# Plot 1  ------
# -----------------
p1 <- ggplot() +
  geom_point(
    data = plot_df %>% filter(),
    aes(x = pos, y = `z`),
    size = 8, alpha = 0.1
  ) +
  facet_grid(
    study + region ~ .,
    labeller = labeller(.rows = custom_labeller),
    scale = "free_y"
  ) +
  geom_point(
    data = plot_df %>% filter(CS1),
    aes(x = pos, y = `z`),
    color = "steelblue",
    size = 8, alpha = 1
  ) +
  xlim(view_win) +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = text_size),
    axis.title = element_text(size = text_size)
  ) +
  xlab("") +
  ylab("z")

# -----------------
# Plot 2------
# -----------------
p2 <- ggplot() +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_text(size = text_size),
    axis.title.x = element_text(size = text_size)
  ) +
  xlim(view_win) +
  ylab("Estimated effect") +
  geom_line(
    data = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5),
    aes_string(y = "fun_plot", x = "x", col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(
    aes(x = view_win[2] + 1000, 
        xend = (haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5))$start[[1]],
        y = 0, yend = 0, col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_segment(
    aes(xend = view_win[1] - 1000,
        x = (haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5))$end[[1]],
        y = 0, yend = 0, col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(
    arrow = arrow(length = unit(0.5, "cm")), 
    aes(x = start, xend = end, y = -0.12 - strand / 10 - 0.02, 
        yend = -0.12 - strand / 10 - 0.02),
    size = 1,
    data = gene_info %>%
      filter(gene_id %in% gene) %>%
      mutate(study = "gene_plot")
  ) +
  geom_text(
    aes(x = (start + end) / 2, 
        y = -0.12 - strand / 10,
        label = gene_name),
    size = 10, 
    data = gene_info %>%
      filter(gene_id %in% gene) %>%
      mutate(study = "gene_plot")
  ) +
  xlab("Phenotype Position") +
  ylab("Estimated\neffect") +
  geom_point(data = sumstat, aes(x = pos - start_distance, y = beta), color = "steelblue", size = 2)

# -----------------
# Plot 3 -----
# -----------------
p3 <- ggplot() +
  facet_grid(
    study ~ .,
    scale = "free_y",
    labeller = labeller(study = c(ROSMAP_DLPFC_haQTL = "car-QTL"))
  ) +
  xlim(view_win) +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = text_size),
    axis.title = element_text(size = text_size)
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, 1), n = 3),
    limits = c(0, 1)
  ) +
  xlab("") +
  ylab("PIP") +
  geom_point(
    data = pip_df %>%
      filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95 == 5),
    aes(x = pos, y = pip, color = as.character(cs_coverage_0.95)),
    alpha = 1, size = 8, color = "steelblue"
  ) +
  xlab("Genotype Position") +
  ylab("PIP") 

# Combine Plots
e <- cowplot::plot_grid(plotlist = list(p1, p3, p2),
                        ncol = 1,
                        align = "v",
                        axis = "tlbr",
                        label_size = 45,  
                        label_fontface = "bold",
                        rel_heights = c(4, 2, 4)
) 

e 
### My code start here ----
plot_df <- data$e[[1]]
haQTL_df <- data$e[[2]]
gene_info <- data$e[[3]]
sumstat <- data$e[[4]]
pip_df <- data$e[[5]]






chr=1
study="AD_Bellenguez_2022"
data_track = plot_df [ which(   plot_df $study == study   ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1   ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black



otAD <- OverlayTrack(trackList=list(    t1,  t2 ))



plotTracks( otAD )




#### CR1 panel ---- 


gene_name="CR1"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black


otCR1 <- OverlayTrack(trackList=list(    t1, t2 ))


plotTracks( otCR1 )




#### CR2 panel ---- 


gene_name="CR2"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black


otCR2 <- OverlayTrack(trackList=list(    t1, t2 ))


plotTracks( otCR2 )



list_track=  list( otAD,
                   otCR1,
                   otCR2  
)

view_win <- c(207317782, 207895513)
plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])





## gene track plot ----
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Initialize AnnotationHub
ah <- AnnotationHub()

# Load the TxDb for GRCh38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

chr=paste0("chr",1)
# Extract the relevant genes and exons in the specified region
region_genes <- genes(txdb,columns = c("tx_id","gene_id"))

# Subset the genes and exons to the region of interest.
region_genes <- subsetByOverlaps(region_genes,
                                 GRanges(seqnames = chr,
                                         ranges = IRanges(view_win[1],
                                                          view_win[2])))
cex <- 0.6
# Create a gene region track for the specified region
gene_track <- GeneRegionTrack(txdb,genome = "hg38",chromosome = chr,
                              pos0 = view_win[1],pos1 =view_win[2],name = "",
                              showId = TRUE,geneSymbol = TRUE,
                              col.axis = "black",col.title = "black",
                              transcriptAnnotation = "symbol",
                              rotation.title = 0,cex.title = cex,
                              col = "salmon",fill = "salmon",
                              background.title = "white")
# Map gene IDs to gene symbols
gene_ids <- unique(unlist(region_genes$gene_id))  # Get unique gene IDs

# Map to gene symbols using org.Hs.eg.db
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID")

  

if(nrow(gene_symbols)>0){
  for(i in 1:length(gene_symbols$ENTREZID)){
    gene_track@range@elementMetadata@listData$id[gene_track@range@elementMetadata@listData$gene == gene_symbols$ENTREZID[i]] <- gene_symbols$SYMBOL[i]
    gene_track@range@elementMetadata@listData$symbol[gene_track@range@elementMetadata@listData$gene == gene_symbols$ENTREZID[i]] <- gene_symbols$SYMBOL[i]
    
    
  }
}

plotTracks(gene_track,
           from =view_win[1],
           to=view_win[2])




res <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_mQTL.chr1_205117782_208795513.fsusie_mixture_normal_top_pc_weights.rds")
res
fsusie_obj_me = res$`chr1:205117782-208795513`$ROSMAP_DLPFC_mQTL$fsusie_result
rm(res)