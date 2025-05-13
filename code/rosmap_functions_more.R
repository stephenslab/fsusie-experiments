# This is a function used to extract the gene annotations from the GTF
# (“gene transfer format”) file. Only the annotated gene transcripts
# for protein-coding genes as defined in the Ensembl/Havana database
# are kept.
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

# Helper function for loading the enrichment results. The "n" argument
# specifies the number of "meta data" columns.  Columns after that are
# treated as the enrichment results. These columns contain only binary
# data (0 or 1) indicating whether or not the genomic feature (genetic
# variant or molecular trait location) is assigned that specific
# annotation.
read_enrichment_results <- function (filename, n) {
  out <- fread(filename,sep = "\t",stringsAsFactors = FALSE,header = TRUE)
  class(out) <- "data.frame"
  out <- transform(out,chr = factor(chr))
  if (ncol(out) > n) {
    cols <- seq(n + 1,ncol(out))
    for (i in cols)
      out[[i]] <- factor(out[[i]])
  }
  return(out)
}

# his function removes CSs so that no two CSs share the same SNP.
create_cs_maps <- function (info_df) {
  info_df %>% 
    arrange(cs) %>%
    group_by(variant_id) %>%
    mutate(cs_set = list(sort(unique(cs))),
           cs_set = sapply(cs_set,paste,collapse = ",")) %>%
    ungroup %>%
    count(cs,cs_set) %>%
    group_by(cs) %>%
    summarize(
	  cs_set = paste(sort(unique(unlist(strsplit(cs_set,",")))),
	                 collapse = ","),.groups = "drop") %>%
          merge_cs_sets
}

# Function written by Hao Sun.
merge_cs_sets <- function (df) {
    
  # df is expected to have at least columns: cs, cs_set
  # cs_set is a comma-separated list of elements belonging to that CS.
  
  # 1. Extract sets: split cs_set by comma
  sets <- strsplit(df$cs_set, ",")
  names(sets) <- df$cs
  
  # 2. Build an element-to-sets map to identify which sets share elements
  element_map <- new.env(hash = TRUE, parent = emptyenv())
  for (s in names(sets)) {
    for (elem in sets[[s]]) {
      if (!exists(elem, envir = element_map)) {
        assign(elem, s, envir = element_map)
      } else {
        assign(elem, c(get(elem, envir = element_map), s), envir = element_map)
      }
    }
  }
  
  # Union-Find (Disjoint Set) Setup
  cs_names <- names(sets)
  parent <- setNames(cs_names, cs_names)
  rank <- setNames(rep(0, length(cs_names)), cs_names)
  
  find_set <- function(x) {
    if (parent[[x]] != x) {
      parent[[x]] <<- find_set(parent[[x]])
    }
    parent[[x]]
  }
  
  union_set <- function(x, y) {
    rx <- find_set(x)
    ry <- find_set(y)
    
    if (rx != ry) {
      if (rank[[rx]] < rank[[ry]]) {
        parent[[rx]] <<- ry
      } else if (rank[[rx]] > rank[[ry]]) {
        parent[[ry]] <<- rx
      } else {
        parent[[ry]] <<- rx
        rank[[rx]] <<- rank[[rx]] + 1
      }
    }
  }
  
  # 3. Union sets that share elements
  # For each element, union all sets that contain it.
  for (elem in ls(element_map)) {
    sets_for_elem <- get(elem, envir = element_map)
    if (length(sets_for_elem) > 1) {
      base_set <- sets_for_elem[1]
      for (other_set in sets_for_elem[-1]) {
        union_set(base_set, other_set)
      }
    }
  }
  
  # 4. Map each cs to its root
  cs_to_root <- data.frame(
    cs = cs_names,
    root = sapply(cs_names, find_set),
    stringsAsFactors = FALSE
  )
  
  # Group by root to get merged sets
  # For each merged component, combine all elements from its sets
  merged_sets <- cs_to_root %>%
    group_by(root) %>%
    summarize(
      merged_sets = list(unique(unlist(sets[cs]))),
      .groups = "drop"
    ) %>%
    mutate(merged_sets_str = sapply(merged_sets, function(x) paste(sort(x), collapse = ",")))%>%select(-merged_sets)
 
  # Return both the mapping of each CS to its root and the merged sets
return(left_join(cs_to_root, merged_sets))
}

# This function is used below to get the sizes of the TADs (in Mb).
get_tad_sizes <- function (tads) {
  tads <- strsplit(tads,"_",fixed = TRUE)
  pos0 <- as.numeric(sapply(tads,"[[",2))
  pos1 <- as.numeric(sapply(tads,"[[",3))
  return((pos1 - pos0)/1e6)
}

# This function is used to summarize the number of CSs per TAD.
get_cs_vs_tad_size <- function (dat) {
  tads <- levels(dat$region)
  out <- data.frame(tad      = tads,
                    tad_size = get_tad_sizes(tads),
                    num_cs   = tapply(dat$cs,dat$region,
                                      function (x) length(unique(x))))
  rownames(out) <- NULL
  return(out)
}

# This function is used to extract the top SNP per location (e.g.,
# CpG) from the association tests.
get_top_snp_per_location <- function (dat) {
  x <- factor(dat$molecular_trait_id)
  qval <- dat$qvalue
  names(qval) <- with(dat,paste(chr,pos,sep = "_"))
  res <- tapply(qval,x,function (x) names(which.min(x)))
  res <- strsplit(res,"_",fixed = TRUE)
  out <- data.frame(chr = factor(sapply(res,"[[",1)),
                    pos = as.numeric(sapply(res,"[[",2)))
  out <- transform(out,chr = factor(chr))
  rownames(out) <- levels(x)
  return(out)
}

# This function adds a column to the SNP results containing the minimum
# distance to the nearest TSS.
add_min_dist_to_tss <- function (dat, genes) {
  n <- nrow(dat)
  dat$min_dist_to_tss <- rep(Inf,n)
  n <- nrow(genes)
  for (i in 1:n) {
    rows <- which(as.character(genes[i,"chromosome"]) == as.character(dat$chr))
	if (length(rows) > 0) {
      if (genes[i,"strand"] == "+")
        d <- genes[i,"start"] - dat[rows,"pos"]
      else
        d <- dat[rows,"pos"] - genes[i,"end"]
	  i <- which(abs(d) < abs(dat[rows,"min_dist_to_tss"]))
  	  if (length(i) > 0) {
	    d    <- d[i]
	    rows <- rows[i]
	    dat[rows,"min_dist_to_tss"] <- d
      }
    }
  }
  return(dat)
}

