# Load libraries
suppressMessages(library(ape))

consensus_blopt_directory <- args[1]
alerax_directory <- args[2]
output_directory <- args[3]
OG_id <- args[4]

final_tree_file <- paste0(consensus_blopt_directory, "/", OG_id, '_consensus50_blopt.treefile')
alerax_consensus_file <- paste0(alerax_directory, "/", OG_id, "/reconciliations/summaries/family_1_consensus_50.newick")
final_tree_file_basename <- gsub("\\..*", "", gsub(".*\\/", "", final_tree_file))
out_tree_file <- paste0(output_directory, '/', final_tree_file_basename, '_supports.treefile')

if (file.exists(out_tree_file)) {
  quit(save = "no")
}

normalize_labels <- function(x) trimws(as.character(x))

# Get descendant tip labels for a node
get_tips <- function(tr, node) {
  if ("Descendants" %in% ls("package:ape")) {
    idx <- ape::Descendants(tr, node, type = "tips")[[1]]
    tr$tip.label[idx]
  } else {
    kids <- function(n) tr$edge[tr$edge[,1] == n, 2]
    collect <- function(n) {
      k <- kids(n)
      if (length(k) == 0) return(n)
      unlist(lapply(k, collect))
    }
    idx <- collect(node)
    tr$tip.label[idx[idx <= length(tr$tip.label)]]
  }
}

# Canonical unrooted bipartition signature (order-independent)
node_signature <- function(tr, node) {
  tips_all <- tr$tip.label
  desc <- sort(get_tips(tr, node))
  comp <- sort(setdiff(tips_all, desc))
  key1 <- paste0(desc, collapse = ";")
  key2 <- paste0(comp, collapse = ";")
  if (key1 <= key2) key1 else key2
}

# Map signature to numeric CCP support
ccp_map_from_tree <- function(tr) {
  if (is.null(tr$node.label)) tr$node.label <- rep(NA_character_, tr$Nnode)
  vals <- suppressWarnings(as.numeric(tr$node.label))
  names(vals) <- sapply((Ntip(tr)+1):(Ntip(tr)+tr$Nnode),
                        function(n) node_signature(tr, n))
  vals
}

lookup_supports <- function(target_tr, lookup_map) {
  sigs <- sapply((Ntip(target_tr)+1):(Ntip(target_tr)+target_tr$Nnode),
                 function(n) node_signature(target_tr, n))
  vals <- unname(lookup_map[sigs])
  as.numeric(vals)
}

# Read trees
final_tr  <- read.tree(final_tree_file)
alerax_tr <- read.tree(alerax_consensus_file)

# Normalize tip labels (trim spaces)
final_tr$tip.label  <- normalize_labels(final_tr$tip.label)
alerax_tr$tip.label <- normalize_labels(alerax_tr$tip.label)

# Align leaves
final_tips <- sort(final_tr$tip.label)
if (!all(final_tips %in% alerax_tr$tip.label)) {
  missing_in_alerax <- setdiff(final_tips, alerax_tr$tip.label)
  stop(sprintf(
    "Tips in the final tree but missing in the AleRax consensus: %s",
    paste(missing_in_alerax, collapse = ", ")
  ))
}
if (!setequal(final_tips, alerax_tr$tip.label)) {
  to_drop <- setdiff(alerax_tr$tip.label, final_tips)
  if (length(to_drop)) {
    message("Pruning extra tips from AleRax tree: ", paste(to_drop, collapse = ", "))
    alerax_tr <- drop.tip(alerax_tr, to_drop)
  }
}

# Map CCP supports
ccp_map  <- ccp_map_from_tree(alerax_tr)
ccp_vals <- lookup_supports(final_tr, ccp_map)

# Write CCP supports to nodes
final_tr$node.label <- ifelse(is.na(ccp_vals), "", as.character(ccp_vals))
write.tree(final_tr, file = out_tree_file)


