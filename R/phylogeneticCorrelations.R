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
#' @param phy an object of class \code{\link[ape]{phylo}}.
#' It must be either a tree with tips having the same names as the columns of \code{object} (including replicates),
#' or a tree such that tip labels match with species names in `col_species`.
#' @param col_species a character vector with same length as columns in the expression matrix,
#' specifying the species for the corresponding column. If left `NULL`, an automatic parsing of species names with sample ids is attempted.
#' @param model the phylogenetic model used to correct for the phylogeny.
#' Must be one of "BM", "lambda", "OUfixedRoot", "OUrandomRoot" or "delta".
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param measurement_error a logical value indicating whether there is measurement error.
#' Default to \code{TRUE}.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param trim the fraction of observations to be trimmed from each end when computing the trimmed mean. Default to 0.15, as in \code{\link[limma]{duplicateCorrelation}}.
#' @param weights a named vector or matrix with weights to be applied on the measurement error term.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param REML Use REML (default) or ML for estimating the parameters.
#' @param ddf_method the method for the computation of the degrees of freedom of the t statistics (before moderation).
#' Default to \code{ddf_method="Satterthwaite"}.
#' If \code{ddf_method="Species"}, then the number of species is taken for the
#' computation of the degrees of freedom,
#' while if \code{ddf_method="Samples"} the total number of individuals is used.
#' @param ... further parameters to be passed
#' to \code{\link[limma]{lmFit}}.
#'
#' @return An object of class \code{TransTree-class},
#' with list components:
#' \itemize{
#' \item \code{tree} the transformed consensus tree
#' \item \code{params} the associated consensus parameters.
#' }
#' This consensus tree defines a correlation structure, and can be passed on to \code{\link[limma]{lmFit}}.
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[phylolm]{phylolm}}, \code{\link[limma]{duplicateCorrelation}}
#'
#' @importFrom methods new
#'
#' @export
#'
#' @importFrom graphics lines title
#' @importFrom stats approxfun lowess model.matrix
#'
phylogeneticCorrelations <- function(object, design = NULL, phy, col_species = NULL,
                                     model = c("BM", "lambda", "OUfixedRoot", "OUrandomRoot", "delta"),
                                     measurement_error = TRUE,
                                     trim = 0.15, weights = NULL, REML = TRUE,
                                     ddf_method = c("Satterthwaite", "Species", "Samples"),
                                     ...) {

  ##################################################################################################
  ## Checks

  ## Expression Matrix
  if (!is.matrix(object)) stop("'object' must be a matrix.")
  y <- limma::getEAWP(object)
  if (!nrow(y$exprs)) stop("expression matrix has zero rows")

  #	Check weights
  if(!is.null(weights)) {
    stop("weights are not allowed with the phylogenetic regression (yet).")
    # message("'weights' will be used in the independent errors.")
    # weights <- limma::asMatrixWeights(weights, dim(y))
    # weights[weights <= 0] <- NA
    # y[!is.finite(weights)] <- NA
  }

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
  if(!is.null(ne)) stop("Coefficients not estimable: ", paste(ne, collapse = " "), "\n")

  ## phylo model
  model <- match.arg(model)

  ## tree
  if (!inherits(phy, "phylo")) stop("object 'phy' must be of class 'phylo'.")
  if (length(phy$tip.label) == ncol(y$exprs)) {
    tree_rep <- phy
  } else {
    if (is.null(col_species)) col_species <- parse_species(phy, colnames(y$exprs))
    tt <- data.frame(species = col_species,
                     id = colnames(y$exprs))
    tree_rep <- addReplicatesOnTree(phy, tt)
    tree_norep <- phy
    phy <- tree_rep
  }
  y_data <- checkParamMatrix(y$exprs, "expression matrix", phy)
  design <- checkParamMatrix(design, "design matrix", phy, transpose = TRUE)

  ## ddf
  ddf_method <- match.arg(ddf_method)

  ##################################################################################################

  tree_model <- get_consensus_tree(y_data, design, phy, model, measurement_error, weights, trim, REML, ddf_method, ...)

  return(tree_model)
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
get_consensus_tree <- function(y_data, design, phy, model, measurement_error, weights, trim, REML, ddf_method, ...) {
  if(!is.null(weights)) stop("weights are not allowed with the phylogenetic regression (yet).")

  if (model == "BM" && !measurement_error) # no parameter to estimate
    return(list(tree = phy,
                params = list(model = "BM",
                              measurement_error = FALSE),
                ddf = rep(nrow(design) - ncol(design), nrow(y_data))))

  # flag_BM_error <- FALSE
  # if (model == "BM" && measurement_error) {
  #   model <- "lambda"
  #   measurement_error <- FALSE
  #   flag_BM_error <- TRUE
  # }

  all_fits <- list()
  for (i in 1:nrow(y_data)) {
    y <- y_data[i, ]
    # w <- weights[i, ]

    data_phylolm <- as.data.frame(cbind(y, design))
    colnames(data_phylolm)[1] <- "expr"
    alpha_bounds <- getBoundsSelectionStrength(phy, 0.0001, 100)
    min_error <- getMinError(phy)
    lower_bounds <- get_lower_bounds(alpha_bounds, min_error, ...)
    upper_bounds <- get_upper_bounds(alpha_bounds, min_error, ...)
    starting_values <- get_starting_values(alpha_bounds, ...)
    dots_args <- get_dots_args(...)

    nafun <- function(e) {
      if (grepl("distinguishable", e)) stop(e)
      return(list(optpar = NA,
                  sigma2 = NA,
                  sigma2_error = NA))
    }
    tmp_fun <- function(...) {
      return(tryCatch(withCallingHandlers(phylolm::phylolm(expr ~ -1 + .,
                                                           data = data_phylolm, phy = phy, model = model,
                                                           measurement_error = measurement_error,
                                                           lower.bound = lower_bounds,
                                                           upper.bound = upper_bounds,
                                                           starting.value = starting_values,
                                                           REML = REML),
                                          warning = function(cond) {
                                            if (grepl(pattern="upper/lower", x = conditionMessage(cond)) && "warning" %in% class(cond)) invokeRestart("muffleWarning")
                                          }),
                      error = nafun))
      # error_weight = weights, ...))
    }
    all_fits[[i]] <- do.call(tmp_fun, dots_args)

  }

  dot_args <- dots(...)
  if (!"medianOU" %in% names(dot_args)) {
    medianOU <- FALSE
  } else {
    medianOU <- dot_args$medianOU
  }

  params <- switch(model,
                   BM = get_consensus_tree_BM(phy, all_fits, measurement_error, trim),
                   lambda = get_consensus_tree_lambda(phy, all_fits, measurement_error, trim),
                   OUfixedRoot = get_consensus_tree_OUfixedRoot(phy, all_fits, measurement_error, trim, medianOU, alpha_bounds),
                   OUrandomRoot = get_consensus_tree_OUrandomRoot(phy, all_fits, measurement_error, trim, medianOU),
                   delta = get_consensus_tree_delta(phy, all_fits, measurement_error, trim))
  params$ddf <- sapply(all_fits, get_ddf(ddf_method), phylo = phy)
  # if (flag_BM_error) {
  #   params$model <- "BM"
  #   params$measurement_error <- TRUE
  # }
  return(params)
}

#' @title Get lambda transformed tree
#'
#' @description
#' Compute the consensus parameters.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
get_consensus_tree_lambda <- function(phy, all_phyfit, measurement_error, trim) {
  if (measurement_error) stop("Measurement error is not allowed with lambda model.")

  all_lambdas <- sapply(all_phyfit, function(x) x$optpar)
  all_lambdas_transform <- atanh(pmax(-1, all_lambdas))
  lambda_mean <- tanh(mean(all_lambdas_transform, trim = trim, na.rm = TRUE))
  tree_model <- phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(tree = tree_model,
              params = list(model = "lambda",
                            measurement_error = measurement_error,
                            lambda = lambda_mean,
                            atanh_lambda = all_lambdas_transform)))
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
get_consensus_tree_BM <- function(phy, all_phyfit, measurement_error, trim) {

  if (!measurement_error) stop("Measurement error should be true here.")

  h_tree <- tree_height(phy)
  all_lambda_error <- sapply(all_phyfit, function(phyfit) get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, h_tree))

  all_lambdas_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_mean <- tanh(mean(all_lambdas_transform, trim = trim, na.rm = TRUE))

  tree_model <- phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(tree = tree_model,
              params = list(model = "BM",
                            measurement_error = measurement_error,
                            lambda_error = lambda_mean,
                            atanh_lambda_error = all_lambdas_transform)))
}

#' @title Get OU transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#' @param median if TRUE, the \code{pracma::geo_median} function is used to take
#' the geometric median.
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
get_consensus_tree_OUfixedRoot <- function(phy, all_phyfit, measurement_error, trim, median = FALSE, alpha_bounds) {

  all_alphas <- sapply(all_phyfit, function(x) x$optpar)
  all_alphas_transform <- log(all_alphas)
  # alpha_mean <- exp(mean(all_alphas_transform, trim = trim, na.rm = TRUE))
  is_min_alpha <- sapply(all_alphas, function(xx) isTRUE(all.equal(xx, alpha_bounds[1], tolerance = (.Machine$double.eps)^(1/3))))
  is_max_alpha <- sapply(all_alphas, function(xx) isTRUE(all.equal(xx, alpha_bounds[2], tolerance = (.Machine$double.eps)^(1/3))))
  non_min_max <- rep(TRUE, length(is_max_alpha)) #!is_min_alpha # & !is_max_alpha
  tree_ind <- rep("treecons", length(non_min_max))
  tree_ind[is_min_alpha] <- "treemin"
  # tree_ind[is_max_alpha] <- "treemax"

  if (!measurement_error) {

    alpha_mean <- exp(mean(all_alphas_transform[non_min_max], trim = trim, na.rm = TRUE))
    return(list(
      tree = list(
        treecons = rescale_tree(phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_mean))$tree),
        treemin = rescale_tree(phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha =  alpha_bounds[1]))$tree),
        treemax = rescale_tree(phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha =  alpha_bounds[2]))$tree)
      ),
      params = list(model = "OUfixedRoot",
                    measurement_error = measurement_error,
                    alpha = alpha_mean,
                    alpha_min = alpha_bounds[1],
                    alpha_max = alpha_bounds[2],
                    log_alpha = all_alphas_transform,
                    tree_ind = tree_ind))
    )
  }

  get_lambda_error_OU <- function(phyfit) {
    tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
                                                 parameters = list(alpha = phyfit$optpar))$tree
    tilde_t <- tree_height(tree_model) / (2 * phyfit$optpar)
    lambda_ou_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
    return(lambda_ou_error)
  }

  ## Use _OU or _OU_cons ? -> cons does not make sense : sigma2 / 2 alpha* t(alpha) is better estimated
  # get_lambda_error_OU_cons <- function(phyfit) {
  #   tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
  #                                                parameters = list(alpha = alpha_mean))$tree
  #   tilde_t <- tree_height(tree_model) / (2 * alpha_mean)
  #   lambda_ou_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
  #   return(lambda_ou_error)
  # }

  ## consensus lambda error
  all_lambda_error <- sapply(all_phyfit, get_lambda_error_OU)
  all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))

  if (!median) {
    alpha_mean <- exp(mean(all_alphas_transform[non_min_max], trim = trim, na.rm = TRUE))
    lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))
  } else {
    ## Geometric median
    all_pars_OU_transform <- cbind(all_alphas_transform, all_lambda_error_transform)
    gmed <- pracma::geo_median(all_pars_OU_transform)
    gmed <- unname(gmed$p)
    alpha_mean <- exp(gmed[1])
    lambda_error_mean <- tanh(gmed[2])
  }

  ## transform tree
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_mean))$tree
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree
  tree_model <- rescale_tree(tree_model)

  # tree_min <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_bounds[1] * 1/10))$tree
  tree_min <- phy # keep BM tree
  tree_min <- phylolm::transf.branch.lengths(tree_min, "lambda", parameters = list(lambda = lambda_error_mean))$tree
  tree_min <- rescale_tree(tree_min)

  tree_max <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_bounds[2]))$tree
  tree_max <- phylolm::transf.branch.lengths(tree_max, "lambda", parameters = list(lambda = lambda_error_mean))$tree
  tree_max <- rescale_tree(tree_max)

  return(list(
    tree = list(
      treecons = tree_model,
      treemin = tree_min,
      treemax = tree_max
    ),
    params = list(model = "OUfixedRoot",
                  measurement_error = measurement_error,
                  alpha = alpha_mean,
                  alpha_min = alpha_bounds[1],
                  alpha_max = alpha_bounds[2],
                  log_alpha = all_alphas_transform,
                  lambda_error = lambda_error_mean,
                  atanh_lambda_error = all_lambda_error_transform,
                  tree_ind = tree_ind))
  )
}

#' @title Get OU transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#' @param median if TRUE, the \code{pracma::geo_median} function is used to take
#' the geometric median.
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
get_consensus_tree_OUrandomRoot <- function(phy, all_phyfit, measurement_error, trim, median = FALSE) {

  all_alphas <- sapply(all_phyfit, function(x) x$optpar)
  all_alphas_transform <- log(all_alphas)
  # alpha_mean <- exp(mean(all_alphas_transform, trim = trim, na.rm = TRUE))

  if (!measurement_error) {

    alpha_mean <- exp(mean(all_alphas_transform, trim = trim, na.rm = TRUE))
    return(list(tree = phylolm::transf.branch.lengths(phy, "OUrandomRoot", parameters = list(alpha = alpha_mean))$tree,
                params = list(model = "OUrandomRoot",
                              measurement_error = measurement_error,
                              alpha = alpha_mean,
                              log_alpha = all_alphas_transform)))
  }

  get_lambda_error_OU <- function(phyfit) {
    tree_model <- phylolm::transf.branch.lengths(phy, "OUrandomRoot",
                                                 parameters = list(alpha = phyfit$optpar))$tree
    tilde_t <- tree_height(tree_model) / (2 * phyfit$optpar)
    lambda_ou_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
    return(lambda_ou_error)
  }

  ## consensus lambda error
  all_lambda_error <- sapply(all_phyfit, get_lambda_error_OU)
  all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))

  if (!median) {
    alpha_mean <- exp(mean(all_alphas_transform, trim = trim, na.rm = TRUE))
    lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))
  } else {
    ## Geometric median
    all_pars_OU_transform <- cbind(all_alphas_transform, all_lambda_error_transform)
    gmed <- pracma::geo_median(all_pars_OU_transform)
    gmed <- unname(gmed$p)
    alpha_mean <- exp(gmed[1])
    lambda_error_mean <- tanh(gmed[2])
  }

  ## transform tree
  tree_model <- phylolm::transf.branch.lengths(phy, "OUrandomRoot", parameters = list(alpha = alpha_mean))$tree
  tree_model$root.edge <- 0
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(tree = tree_model,
              params = list(model = "OUrandomRoot",
                            measurement_error = measurement_error,
                            alpha = alpha_mean,
                            lambda_error = lambda_error_mean,
                            log_alpha = all_alphas_transform,
                            atanh_lambda_error = all_lambda_error_transform)))
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
get_consensus_tree_delta <- function(phy, all_phyfit, measurement_error, trim) {

  all_deltas <- sapply(all_phyfit, function(x) x$optpar)
  all_deltas_transform <- log(all_deltas)
  delta_mean <- exp(mean(all_deltas_transform, trim = trim, na.rm = TRUE))

  if (!measurement_error) {

    return(list(tree = phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = delta_mean))$tree,
                params = list(model = "delta",
                              measurement_error = measurement_error,
                              delta = delta_mean,
                              log_delta = all_deltas_transform)))
  }

  get_lambda_error_delta <- function(phyfit) {
    tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = phyfit$optpar))$tree
    tilde_t <- tree_height(tree_model)
    lambda_delta_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
    return(lambda_delta_error)
  }

  all_lambda_error <- sapply(all_phyfit, get_lambda_error_delta)
  all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))

  ## transform tree
  tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = delta_mean))$tree
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(tree = tree_model,
              params = list(model = "delta",
                            measurement_error = measurement_error,
                            delta = delta_mean,
                            lambda_error = lambda_error_mean,
                            log_delta = all_deltas_transform,
                            atanh_lambda_error = all_lambda_error_transform)))
}

check.consensus_tree <- function(consensus_tree, model, measurement_error) {
  if (consensus_tree$params$model != model) {
    stop(paste0("The consensus tree was computed with the ", consensus_tree$params$model, " model, but the model used is ", model, ". Please check which model is the correct one."))
  }
  if (consensus_tree$params$measurement_error != measurement_error) {
    stop(paste0("The consensus tree was computed with `measurement_error=", consensus_tree$params$measurement_error, "`, but you specified `measurement_error=", measurement_error, "`. Please check which model is the correct one."))
  }
}
