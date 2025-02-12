rm(list=ls())

library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
path= getwd()
data = readRDS(paste0(  path,"/data/fig_4_data/Fig4_data.rds") )
extract_snp_position <- function(snp_string) {
  # Split the input string by ':'
  parts <- unlist(strsplit(snp_string, ":"))
  
  # Extract the position
  position <- as.numeric(parts[2])
  
  return(position)
}
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



##" cas 1 ---- 
plot_df <- data$d[[1]]
haQTL_df <- data$d[[2]]
mQTL_df <- data$d[[3]]
gene_info <- data$d[[4]]
sumstat <- data$d[[5]]
pip_df <- data$d[[6]]


# Parameters
view_win <- c(4759843, 5000000)
text_size <- 20
gene <- c("ENSG00000029725", "ENSG00000161929", "ENSG00000108556")
custom_labeller <- function(x) {
  x %>% 
    gsub("DeJager_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", "\n", .)             
}
# Function to extract SNP position from a given notation




plot_df <- data$f[[1]]
haQTL_df <- data$f[[2]]
MSBB_df <- data$f[[3]]
gene_info <- data$f[[4]]
sumstat <- data$f[[5]]
pip_df <- data$f[[6]]
QTL_data <- data$f[[7]] 
# Custom Labeller
custom_labeller <- function(x) {
  x %>% 
    gsub("mega_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", " \n ", .)             
}

# Parameters
text_size <- 20
view_win <- c(5.12e7, 5.16e7)
gene <- c("ENSG00000139629", "ENSG00000050438")

# -----------------
# Plot 1
# -----------------
p1 <- ggplot() +
  geom_point(
    data = plot_df %>% filter(),
    aes(x = pos, y = `z`),
    size = 8, alpha = 0.1
  ) +
  facet_grid(
    study + region ~ .,
    labeller = labeller(.rows =custom_labeller),
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
  ylab("Z")

# -----------------
# Plot 2
# -----------------
p2 <- ggplot() +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_text(size = text_size),    
    axis.text.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size)
  ) +
  xlim(view_win) +
  geom_line(
    data = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5),
    aes_string(y = "fun_plot", x = "x", col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_line(
    data = MSBB_df %>% mutate(study = "dmr QTL effect") %>% filter(CS == 14),
    aes_string(y = "fun_plot", x = "x", col = "CS"),
    size = 4, col = "maroon"
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
  ylab("Estimated\neffect")  +
  geom_point(data =  sumstat %>% filter(ha), aes(x = pos - start_distance, y = beta), color = "steelblue", size = 2) +
  geom_point(data =  sumstat %>% filter(!ha), aes(x = pos - start_distance, y = beta), color = "maroon", size = 2)

# -----------------
# Plot 3
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
    axis.text.x = element_text(size = text_size),
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
  geom_point(
    data = QTL_data %>%
      filter(variant_id == "chr12:51362485:T:C", study == "MSBB_mQTL") %>%
      mutate(study = "dmr-QTL") %>%
      separate(col = variant_id, into = c("chrom", "pos"), remove = FALSE) %>%
      mutate(pos = as.numeric(pos)),
    aes(x = pos, y = pip, color = as.character(cs_coverage_0.95)),
    alpha = 1, size = 8, color = "maroon"
  ) +
  xlab("Genotype Position") +
  ylab("PIP") 

# Combine Plots
f <- cowplot::plot_grid(plotlist = list(p1, p3, p2),
                        ncol = 1,
                        align = "v",
                        axis = "tlbr",
                        label_fontface = "bold",
                        rel_heights = c(5, 2, 4)
) 

f

### AD GWAS panel -----

chr =  paste("chr", 12, sep = "")

 

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






#### CHRNE panel ---- 


chr =  paste("chr", 17, sep = "")

gene_name="CHRNE"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"
data_track_CS4 =plot_df [ which(  plot_df$study == study & plot_df$CS4& plot_df$region==gene_name ),]# plot_df %>% filter(CS4, study == "DLPFC_DeJager_eQTL") #"steelblue"



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
                type = "p", col = "maroon", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black


t3= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS4$pos , end = data_track_CS4$pos )),
                data = matrix(data_track_CS4$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue",
                cex=1.5,# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black


otCHRNE <- OverlayTrack(trackList=list(    t1, t3 ,t2 ))


plotTracks( otCHRNE )










#### RABEP1 panel ---- 


chr =  paste("chr", 17, sep = "")

gene_name="RABEP1"
study="Mic_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"
data_track_CS4 =plot_df [ which(  plot_df$study == study & plot_df$CS4& plot_df$region==gene_name ),]# plot_df %>% filter(CS4, study == "DLPFC_DeJager_eQTL") #"steelblue"



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
                type = "p", col = "#2E8B57", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black



otRABEP1 <- OverlayTrack(trackList=list(    t1,  t2 ))


plotTracks( otRABEP1 )





