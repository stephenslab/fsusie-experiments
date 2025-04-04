# Written by Hao Sun.
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
