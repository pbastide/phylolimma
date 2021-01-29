#' @title Add replicates to a tree
#'
#' @description
#' Utility function to add replicates to a tree, as tips with zero length branches.
#'
#' @param tree A phylogenetic tree with n tips.
#' @param traits A data frame containing at least two columns,
#' one with sample ids, and on with species names for each samples.
#' @param species Name of the column containing species names. Default to "species".
#' @param id Name of the column containing samples ids. Default to "id".
#' @param eps A small number to add to terminal branch lengths to avoid true zeros. Default to \code{.Machine$double.eps}.
#'
#' @return A phylogenetic tree with as many tips as the number of rows in \code{traits},
#'  and clusters of tips with zero branch lengths corresponding to replicates.
#'
#' @export
#'
addReplicatesOnTree <- function(tree, traits, species = "species", id = "id", eps = .Machine$double.eps) {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is needed for function 'add_replicates'.", call. = FALSE)
  }
  if (is.null(traits[[id]])) {
    stop("The `traits` data frame should contain a column named ", id, " with sample ids. Plead adjust argument `id` accordingly.")
  }
  if (length(traits[[id]]) != length(unique(traits[[id]]))){
    stop("The samples ids in column named ", id, " should be unique identifiers of the samples.")
  }
  if (is.null(traits[[species]])) {
    stop("The `traits` data frame should contain a column named ", species, " with species names for each sample. Plead adjust argument `species` accordingly.")
  }
  data_tree_cor <- match(traits[[species]], tree$tip.label)
  if (anyNA(data_tree_cor)) {
    # Species in data NOT in the tree
    stop("Species '", paste(unique(traits[[species]][is.na(data_tree_cor)]), collapse = "', '"), "' are in the data but not in the tree. Please remove them from the data before proceding (or use a tree that include them)." )
  }
  tree_data_cor <- match(tree$tip.label, traits[[species]])
  if (anyNA(tree_data_cor)) {
    # Species in tree NOT in the data
    warning("Species '", paste(unique(tree$tip.label[is.na(tree_data_cor)]), collapse = "', '"), "' are in the tree but not in the data. They will be droped from the final tree." )
  }
  ## Make tree
  tree_rep <- tree
  # Add replicates
  for (tip_label in tree$tip.label) {
    all_rep_ids <- traits[[id]][traits[[species]] == tip_label]
    for (rep_id in all_rep_ids) {
      tree_rep <- phytools::bind.tip(tree_rep, tip.label = rep_id,
                                     where = which(tree_rep$tip.label == tip_label))
    }
  }
  # Remove original tips
  tree_rep <- ape::drop.tip(tree_rep, tree$tip.label)
  # No true zeros
  tree_rep$edge.length[tree_rep$edge[, 2] %in% 1:length(tree_rep$tip.label)] <- tree_rep$edge.length[tree_rep$edge[, 2] %in% 1:length(tree_rep$tip.label)] + eps
  # result
  return(tree_rep)
}

