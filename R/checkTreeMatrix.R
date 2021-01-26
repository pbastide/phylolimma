#' @title Check Matrix Parameter
#'
#' @description
#' Check that the parameters are compatible with the tree. Throws an error if not.
#'
#' @param x matrix of parameters being tested.
#' @param name name of the parameter.
#' @param tree A phylogenetic tree with n tips.
#' @param transpose Should the transpose of the matrix be taken ? Default to FALSE.
#'
#' @keywords internal
#'
checkParamMatrix <- function(x, name, tree, transpose = FALSE) {
  N <- length(tree$tip.label)

  if (!transpose) {
    if (ncol(x) != N) {
      stop(paste0("`", name, "` should have as many columns as the number of taxa in the tree."))
    }
    if ((is.null(tree$tip.label) || is.null(colnames(x)))){
      stop(paste0("`", name, "` and/or the tips of the phylogeny are not named. I could not check for consistency. Please give names consistent names to the tree tip labels and the column names of matrix `", name, "` to avoid any ambiguity."))
    } else {
      if (!all(tree$tip.label == colnames(x))){
        # Match
        tree_data_cor <- match(tree$tip.label, colnames(x))
        data_tree_cor <- match(colnames(x), tree$tip.label)
        if (anyNA(tree_data_cor)) {
          # Species in the tree NOT in data
          stop("Species '", tree$tip.label[is.na(tree_data_cor)], "' is in the tree but not in the data.")
        }
        if (anyNA(tree_data_cor)) {
          # Species in data NOT in the tree
          stop("Species '", colnames(x)[is.na(data_tree_cor)], "' is in the data but not in the tree.")
        }
        if (length(unique(tree_data_cor)) != length(tree$tip.label)){
          stop(paste0("`", name, "` names do not match the tip labels."))
        }
        warning(paste0("`", name, "` was not sorted in the correct order, when compared with the tips label. I am re-ordering it."))
        x <- x[, tree_data_cor, drop = FALSE]
      }
    }
  } else {
    if (nrow(x) != N) {
      stop(paste0("`", name, "` should have as many rows as the number of taxa in the tree."))
    }
    if ((is.null(tree$tip.label) || is.null(rownames(x)))){
      stop(paste0("`", name, "` and/or the tips of the phylogeny are not named. I could not check for consistency. Please give names consistent names to the tree tip labels and the row names of matrix `", name, "` to avoid any ambiguity."))
    } else {
      if (!all(tree$tip.label == rownames(x))){
        # Match
        tree_data_cor <- match(tree$tip.label, rownames(x))
        data_tree_cor <- match(rownames(x), tree$tip.label)
        if (anyNA(tree_data_cor)) {
          # Species in the tree NOT in data
          stop("Species '", tree$tip.label[is.na(tree_data_cor)], "' is in the tree but not in the design.")
        }
        if (anyNA(tree_data_cor)) {
          # Species in data NOT in the tree
          stop("Species '", rownames(x)[is.na(data_tree_cor)], "' is in the design but not in the tree.")
        }
        if (length(unique(tree_data_cor)) != length(tree$tip.label)){
          stop(paste0("`", name, "` names do not match the tip labels."))
        }
        warning(paste0("`", name, "` was not sorted in the correct order, when compared with the tips label. I am re-ordering it."))
        x <- x[tree_data_cor, , drop = FALSE]
      }
    }
  }
  return(x)
}
