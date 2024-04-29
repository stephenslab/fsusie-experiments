# Extract the gene annotations from the GTF ("gene transfer format")
# file. Here we keep only the annotated gene transcripts for
# protein-coding genes as defined in the Ensembl/Havana database.
get_gene_annotations <- function (gene_file) {
  out <- fread(file = gene_file,sep = "\t",header = FALSE,skip = 1)
  class(out) <- "data.frame"
  names(out) <- c("chromosome","source","feature","start","end","score",
                    "strand","frame","attributes")
  out <- out[c("chromosome","source","feature","start","end","strand",
               "attributes")]
  out <- transform(out,
                   chromosome = factor(chromosome),
                   source     = factor(source),
                   feature    = factor(feature),
                   strand     = factor(strand))
  out <- subset(out,
                source == "ensembl_havana" &
                feature == "transcript")
  out <-
    transform(out,
      ensembl   = sapply(strsplit(attributes,";"),
                         function (x) substr(x[[1]],10,24)),
      gene_type = sapply(strsplit(attributes,";"),
                         function (x) substr(x[[3]],13,nchar(x[[3]]) - 1)),
      gene_name = sapply(strsplit(attributes,";"),
                         function (x) substr(x[[4]],13,nchar(x[[4]]) - 1)))
  out <- out[-7]
  out <- transform(out,gene_type = factor(gene_type))
  out <- subset(out,gene_type == "protein_coding")
  rownames(out) <- NULL
  return(out)
}

# Create a histogram of the region sizes in Megabases (Mb).
region_sizes_histogram <- function (regions, pips, font_size = 10) {
  regions$pos_min <- sapply(pips,function (x) min(x$pos,na.rm = TRUE))
  regions$pos_max <- sapply(pips,function (x) max(x$pos,na.rm = TRUE))
  regions <- transform(regions,size_bp = pos_max - pos_min)
  return(ggplot(regions,aes(size_bp/1e6)) +
         geom_histogram(color = "white",fill = "darkblue",bins = 64) +
         scale_x_continuous(breaks = seq(0,50,5)) +
         labs(x = "size of region (Mb)",
              y = "number of regions") +
         theme_cowplot(font_size = font_size))
}

# Create a histogram of the region sizes in number of SNPs.
num_snps_histogram <- function (regions, font_size = 10) {
  return(ggplot(susie$regions,aes(num_snps)) +
         geom_histogram(color = "white",fill = "darkblue",bins = 64) +
         scale_x_continuous(breaks = seq(0,1e5,1e4)) +
         labs(x = "number of SNPs",
              y = "number of regions") +
         theme_cowplot(font_size = font_size))
}
         
