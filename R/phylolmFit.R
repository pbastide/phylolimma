#' @title Phylogenetic COmparative Method using LIMMA
#'
#' @description
#' This function applies \code{\link[limma]{lmFit}} to the normalized data,
#' in order to take the phylogeny into account.
#' TODO: explain more.
#'
#' @param object 	A matrix data object containing normalized expression values,
#' with rows corresponding to genes and columns to samples (species).
#' @param design the design matrix of the experiment,
#' with rows corresponding to samples and columns to coefficients to be estimated.
#' Defaults to the unit vector (intercept).
#' @param phy an object of class \code{\link[ape]{phylo}},
#' with tips having the same names as the columns of \code{object}.
#' @param model the phylogenetic model used to correct for the phylogeny.
#' Must be one of "BM", "lambda", "OUfixedRoot" or "delta".
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param measurement_error a logical value indicating whether there is measurement error.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param ... further parameters to be passed
#' to \code{\link[limma]{lmFit}} or \code{\link[phylolm]{phylolm}}.
#'
#' @return An object of class \code{\link[limma]{MArrayLM-class}},
#' with list components \code{coefficients}, \code{stdev.unscaled},
#' \code{sigma} and \code{df.residual}.
#' These quantities take the phylogenetic model into account.
#' The object can be passed to \code{\link[limma]{eBayes}}.
#'
#' @importFrom methods new
#'
#' @export
#'
phylolmFit <- function(object, design = NULL, phy,
                       model = c("BM", "lambda", "OUfixedRoot", "delta"),
                       measurement_error = FALSE, ...) {

  ## Check unused parameters
  dot_args <- dots(...)
  if ("ndups" %in% names(dot_args) && dot_args$ndups != 1) stop("'ndups' can only be '1' in 'phylolmFit' (for now).")
  if ("spacing" %in% names(dot_args) && dot_args$spacing != 1) stop("'spacing' can only be '1' in 'phylolmFit' (for now).")
  if ("weights" %in% names(dot_args) && !is.null(dot_args$weights)) stop("'weights' can only be 'null' in 'phylolmFit' (for now).")
  if ("method" %in% names(dot_args) && dot_args$method != "ls") stop("'method' can only be 'ls' in 'phylolmFit' (for now).")
  if ("correlation" %in% names(dot_args)) stop("'correlation' is not used in 'phylolmFit' (for now).")
  if ("block" %in% names(dot_args) && !is.null(dot_args$block)) stop("'block' can only be 'null' in 'phylolmFit' (for now).")

  ## Expression Matrix
  if (!is.matrix(object)) stop("'object' must be a matrix.")
  y <- limma::getEAWP(object)
  if (!nrow(y$exprs)) stop("expression matrix has zero rows")

  ##	Check design matrix
  if(is.null(design)) design <- y$design
  if(is.null(design)) {
    design <- matrix(1, ncol(y$exprs), 1)
    rownames(design) <- phy$tip.label
  } else {
    design <- as.matrix(design)
    if(mode(design) != "numeric") stop("design must be a numeric matrix")
    if(nrow(design) != ncol(y$exprs)) stop("row dimension of design doesn't match column dimension of data object")
  }
  ne <- limma::nonEstimable(design)
  if(!is.null(ne)) cat("Coefficients not estimable:", paste(ne, collapse = " "), "\n")

  ## phylo model
  # if (model != "BM") stop("'modelphy' can only be 'BM' (for now).")
  model <- match.arg(model)

  ## tree
  if (!inherits(phy, "phylo")) stop("object 'phy' must be of class 'phylo'.")
  y_data <- checkParamMatrix(y$exprs, "expression matrix", phy)
  design <- checkParamMatrix(design, "design matrix", phy, transpose = TRUE)

  ## phylo model
  C_tree_params <- get_chol_tree(y_data,  design, phy, model, measurement_error, ...)
  C_tree <- C_tree_params$C_tree

  ## Transform design and data
  design_trans <- transform_design_tree(C_tree, design)
  y_trans <- t(transform_data_tree(C_tree, y_data))

  ## Apply lmFit
  resLmFit <- lmFitLimma(y_trans, design_trans, ...)

  ## Format
  if (model == "BM" && !measurement_error) {
    resFitFormat <- resLmFit
  } else {
    resFitFormat <- new("PhyloMArrayLM",
                        list(coefficients = do.call(rbind, resLmFit["coefficients", ]),
                             sigma = do.call(c, resLmFit["sigma", ]),
                             stdev.unscaled = do.call(rbind, resLmFit["stdev.unscaled", ]),
                             df.residual = do.call(c, resLmFit["df.residual", ])))
  }

  ## Add slots
  resFitFormat$phy <- phy
  resFitFormat$modelphy <- model
  resFitFormat$measurement_error <- measurement_error
  resFitFormat$phy_trans <- C_tree_params$trans_tree
  resFitFormat$optpar <- C_tree_params$optpar
  resFitFormat$lambda_error <- C_tree_params$lambda_error
  resFitFormat$sigma2_phy <- C_tree_params$sigma2_phy
  resFitFormat$sigma2_error <- C_tree_params$sigma2_error

  ## Result
  return(resFitFormat)
}

#' @title Capture dot arguments
#'
#' @description http://adv-r.had.co.nz/Computing-on-the-language.html#capturing-dots
#'
#' @param ... dots arguments to be captured
#'
#' @return a named list of the arguments in ...
#'
#' @keywords internal
#'
dots <- function(...) {
  eval(substitute(alist(...)))
}

#' @title Fit using limma
#'
#' @description Fit using limma
#'
#'
#' @param y_trans A matrix data object containing normalized
#' and phylogeny transformed expression values,
#' with rows corresponding to genes and columns to samples (species).
#' @param design_trans the phylogeny transformed design matrix of the experiment,
#' with rows corresponding to samples and columns to coefficients to be estimated.
#' Defaults to the unit vector (intercept).
#' @param ... further parameters to be passed to \code{\link[limma]{lmFit}}.
#'
#' @return A list with all the results.
#'
#' @keywords internal
#'
lmFitLimma <- function(y_trans, design_trans, ...) {
  if (!is.list(design_trans)) {
    return(limma::lmFit(object = y_trans, design = design_trans, ...))
  }
  argsfun <- c(ndups = 1, spacing = 1, block = NULL, weights = NULL, method = "ls", dots(...))
  all_res_fit <- mapply(limma::lmFit,
                        object = lapply(seq_len(nrow(y_trans)), function(i) y_trans[i,]),
                        design = design_trans,
                        MoreArgs = argsfun)
  all_res_fit
}

#' @title Get Tree Normalizing Inverse Cholesky
#'
#' @description
#' Compute the whitening cholesky matrix.
#'
#' @param y_data 	A matrix data object containing normalized expression values,
#' with rows corresponding to genes and columns to samples (species).
#' @inheritParams phylolmFit
#'
#' @return The (list of) cholesky matrix of the tree structure.
#'
#' @keywords internal
#'
get_chol_tree <- function(y_data, design, phy, model, measurement_error, ...) {
  if (!measurement_error && model == "BM") { ## Easy case, not fit necessary
    C_tree <- ape::vcv(phy)
    C_tree_chol <- chol(C_tree)
    return(list(C_tree = C_tree_chol,
                trans_tree = phy,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = NA,
                sigma2_error = 0))
  } else {
    C_tree_chol_and_params <- apply(y_data, 1,
                                    get_C_tree, design, phy, model, measurement_error, ...)
    C_tree_chol_and_params <- format_list(C_tree_chol_and_params)
    return(C_tree_chol_and_params)
  }
}

format_list <- function(C_tree_chol_and_params) {
  return(list(C_tree = lapply(C_tree_chol_and_params, function(x) x$C_tree),
              trans_tree = lapply(C_tree_chol_and_params, function(x) x$trans_tree),
              optpar = sapply(C_tree_chol_and_params, function(x) x$optpar),
              lambda_error = sapply(C_tree_chol_and_params, function(x) x$lambda_error),
              sigma2_phy = sapply(C_tree_chol_and_params, function(x) x$sigma2_phy),
              sigma2_error = sapply(C_tree_chol_and_params, function(x) x$sigma2_error)))
}

#' @title Get Tree Normalizing Inverse Cholesky
#'
#' @description
#' Compute the whitening cholesky matrix.
#'
#' @param y 	A vector data containing normalized expression values for one gene.
#' @inheritParams phylolmFit
#'
#' @return The cholesky matrix of the tree structure.
#'
#' @keywords internal
#'
get_C_tree <- function(y, design, phy, model, measurement_error, ...) {
  trans_tree_params <- transform_tree_phylolm(y, design, phy, model, measurement_error, ...)
  tree_model <- trans_tree_params$tree_model
  C_tree <- ape::vcv(tree_model)
  C_tree_chol <- chol(C_tree)
  return(list(C_tree = C_tree_chol,
              trans_tree = tree_model,
              optpar = trans_tree_params$optpar,
              lambda_error = trans_tree_params$lambda_error,
              sigma2_phy = trans_tree_params$sigma2_phy,
              sigma2_error = trans_tree_params$sigma2_error))
}

#' @title Get transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_phylolm <- function(y, design, phy, model, measurement_error, ...) {
  if (model == "BM" && !measurement_error) return(phy) # no transformation needed
  data_phylolm <- as.data.frame(cbind(y, design))
  colnames(data_phylolm)[1] <- "expr"
  fplm <- phylolm::phylolm(expr ~ -1 + ., data = data_phylolm, phy = phy, model = model, measurement_error = measurement_error, ...)
  phy_trans_params <- switch(model,
                             BM = transform_tree_model_BM(phy, fplm, measurement_error),
                             lambda = transform_tree_model_lambda(phy, fplm, measurement_error),
                             OUfixedRoot = transform_tree_model_OUfixedRoot(phy, fplm, measurement_error),
                             delta = transform_tree_model_delta(phy, fplm, measurement_error))
  return(phy_trans_params)
}

#' @title Get lambda transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_lambda <- function(phy, phyfit, measurement_error) {
  if (measurement_error) stop("Measurement error is not allowed with lambda model.")
  return(list(tree_model = phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = phyfit$optpar))$tree,
              optpar = NA,
              lambda_error = 1,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = 0))
}

#' @title Get BM transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_BM <- function(phy, phyfit, measurement_error) {
  if (!measurement_error) stop("Measurement error should be true here.")
  lambda_error <- phyfit$sigma2 / (phyfit$sigma2_error / max(ape::vcv(phy)) + phyfit$sigma2)
  return(list(tree_model = phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_error))$tree,
              optpar = NA,
              lambda_error = lambda_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Get OU transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_OUfixedRoot <- function(phy, phyfit, measurement_error) {
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = phyfit$optpar))$tree
  if (!measurement_error) {
    return(list(tree_model = tree_model,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = phyfit$sigma2,
                sigma2_error = 0))
  }
  tilde_t <- max(ape::vcv(tree_model)) / (2 * phyfit$optpar)
  lambda_ou_error <- phyfit$sigma2 * tilde_t / (phyfit$sigma2_error + phyfit$sigma2 * tilde_t)
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_ou_error))$tree
  return(list(tree_model = tree_model,
              optpar = phyfit$optpar,
              lambda_error = lambda_ou_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Get delta transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_delta <- function(phy, phyfit, measurement_error) {
  tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = phyfit$optpar))$tree
  if (!measurement_error) {
    return(list(tree_model = tree_model,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = phyfit$sigma2,
                sigma2_error = 0))
  }
  tilde_t <- max(ape::vcv(tree_model))
  lambda_delta_error <- phyfit$sigma2 * tilde_t / (phyfit$sigma2_error + phyfit$sigma2 * tilde_t)
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_delta_error))$tree
  return(list(tree_model = tree_model,
              optpar = phyfit$optpar,
              lambda_error = lambda_delta_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Transform design matrix
#'
#' @description
#' Multiply by inverse cholesky to whiten the data
#'
#' @param C_tree Cholesky of the tree structure obtained through \code{\link{get_chol_tree}}
#' @param design design matrix
#'
#' @return transformed design matrix
#'
#' @keywords internal
#'
transform_design_tree <- function(C_tree, design) {
  if (!is.list(C_tree)) return(transform_design_one_tree(C_tree, design))
  return(lapply(C_tree, transform_design_one_tree, design))
}

transform_design_one_tree <- function(C_tree, design, transpose = FALSE) {
  return(backsolve(C_tree, design, transpose = TRUE))
}

#' @title Transform data matrix
#'
#' @description
#' Multiply by inverse cholesky to whiten the data
#'
#' @param C_tree Cholesky of the tree structure obtained through \code{\link{get_chol_tree}}
#' @param y_data data matrix
#'
#' @return transformed data matrix
#'
#' @keywords internal
#'
transform_data_tree <- function(C_tree, y_data) {
  if (!is.list(C_tree)) return(transform_design_one_tree(C_tree, t(y_data)))
  return(mapply(transform_design_one_tree,
                C_tree,
                lapply(seq_len(nrow(y_data)), function(i) y_data[i,])))
}
