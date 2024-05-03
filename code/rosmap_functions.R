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

# Extract the base-pair positions from the SNP ids.
get_pos_from_id <- function (ids)
  as.numeric(sapply(strsplit(ids,":"),"[",2))

# Get the TSS and strand (+/-) for each region (gene) from the "genes"
# table. This code is based on the "gtf_to_tss_bed" function from
# https://github.com/broadinstitute/pyqtl/blob/master/qtl/io.py.
add_tss_to_regions <- function (regions, genes) {
  rownames(regions) <- regions$region_name
  regions$tss <- as.numeric(NA)
  regions$strand <- as.character(NA)
  for (i in regions$region_name) {
    j <- which(genes$ensembl == i)
    if (length(j) == 1) {
      regions[i,"strand"] <- as.character(genes[j,"strand"])
      if (genes[j,"strand"] == "+")
        regions[i,"tss"] <- genes[j,"start"]
      else
        regions[i,"tss"] <- genes[j,"end"]
    }
  }
  rownames(regions) <- NULL
  return(transform(regions,strand = factor(strand)))
}

# From a "CS" data frame containing the information about the credible
# sets by region, generate a summary of the CS sizes by region.
get_cs_sizes_by_region <- function (cs) {
  n <- nlevels(cs$region)
  out <- vector("list",n)
  region_names <- levels(cs$region)
  names(out) <- region_names
  for (i in region_names) {
    res <- table(factor(subset(cs,region == i)$cs))
    names(res) <- paste0(i,"-CS",names(res))
    out[[i]] <- res
  }
  return(out)
}

# Compute the distance to the TSS weighted by the PIPs, only for SNPs
# that are in CSs. Note that TSS and strand information needs to be
# added first to "regions" using function add_tss_to_regions.
compute_weighted_distance_to_tss <- function (regions, cs, bins) {
  n <- nrow(regions)
  rownames(regions) <- regions$region_name
  total_counts      <- rep(0,length(bins) - 1)
  
  # Repeat for each region.
  for (i in regions$region_name) {
    tss    <- regions[i,"tss"]
    strand <- regions[i,"strand"]
    if (!is.na(tss)) {
      dat <- subset(cs,region == i)
      d   <- tss - dat$pos
      if (strand == "-") 
        d <- -d
      d <- cut(d,bins)
      counts <- tapply(dat$pip,d,function (x) sum(x,na.rm = TRUE))
      counts[is.na(counts)] <- 0
      total_counts <- total_counts + counts
    }
  }
  
  return(total_counts)
}

# Put together a data frame, extracted from the "pips" data frame,
# containing only the information about the SNPs with larger PIPs
# according to some specified threshold (we'll call these
# "high-confidence SNPs").
get_highconf_snps <- function (pips, level = 0.95) {
  n <- length(pips)
  out <- vector("list",n)
  region_names <- names(pips)
  names(out) <- names(pips)
  for (i in names(pips))
    out[[i]] <- subset(pips[[i]],pip > level)
  out <- do.call(rbind,out)
  rownames(out) <- NULL
  return(out)
}

# Create a histogram of the region sizes in Megabases (Mb).
region_sizes_histogram <- function (regions, pips, max_mb = Inf,
                                    font_size = 10) {
  regions$pos_min <- sapply(pips,function (x) min(x$pos,na.rm = TRUE))
  regions$pos_max <- sapply(pips,function (x) max(x$pos,na.rm = TRUE))
  regions <- transform(regions,size_bp = pos_max - pos_min)
  if (is.finite(max_mb))
    cat(sum(regions$size_bp > 1e6*max_mb),"regions are larger than",
        max_mb,"Mb.\n")
  regions <- subset(regions,size_bp <= 1e6*max_mb)
  return(ggplot(regions,aes(size_bp/1e6)) +
         geom_histogram(color = "white",fill = "darkblue",bins = 64) +
         scale_x_continuous(breaks = seq(0,100,1)) +
         labs(x = "size of region (Mb)",
              y = "number of regions") +
         theme_cowplot(font_size = font_size))
}

# Create a histogram of the region sizes in number of SNPs.
num_snps_histogram <- function (regions, max_snps = Inf, font_size = 10) {
  if (is.finite(max_snps))
    cat(sum(regions$num_snps > max_snps),"regions have more than",max_snps,
        "SNPs.\n")
  regions <- subset(regions,num_snps <= max_snps)
  return(ggplot(regions,aes(num_snps/1000)) +
         geom_histogram(color = "white",fill = "darkblue",bins = 64) +
         scale_x_continuous(breaks = seq(0,100,10)) +
         labs(x = "number of SNPs x 1,000",
              y = "number of regions") +
         theme_cowplot(font_size = font_size))
}

# Create a histogram of the CS sizes in which "cs_sizes" is an output
# from function get_cs_sizes_by_region.
cs_sizes_histogram <- function (cs_sizes, x_breaks = c(1,10,100,1000),
                                font_size = 10) {
  pdat <- data.frame(cs_size = unlist(cs_sizes))
  return(ggplot(pdat,aes(cs_size)) +
         geom_histogram(color = "white",fill = "darkblue",bins = 64) +
         scale_x_continuous(trans = "log10",breaks = x_breaks) +
         labs(x = "number of SNPs",
              y = "number of CSs") +
         theme_cowplot(font_size = font_size))
}
