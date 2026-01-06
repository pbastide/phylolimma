#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
NULL

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
#' Must be one of "BM", "lambda" or "OUfixedRoot".
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param measurement_error a logical value indicating whether there is measurement error.
#' Default to \code{TRUE}.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param trim a vector of size two, with the fraction of observations to be trimmed from the lower and upper ends of `atanh(all.lambdas)` and `atanh(rho)` when computing the trimmed mean. If a single value is provided, it is recycled as a vector of size two. Default to `c(0.25, 0.05)`. See also the `trim` argument in \code{\link[limma]{duplicateCorrelation}}.
#' @param weights a named vector or matrix with weights to be applied on the measurement error term.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param REML Use REML (default) or ML for estimating the parameters.
#' @param ncores number of cores to use for parallel computation. Default to 1 (no parallel computation).
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
#' @importFrom stats approxfun lowess model.matrix uniroot complete.cases
#'
phylogeneticCorrelations <- function(object, design = NULL, phy, col_species = NULL,
                                     model = c("BM", "lambda", "OUfixedRoot"),
                                     measurement_error = TRUE,
                                     trim = c(0.25, 0.05), weights = NULL, REML = TRUE,
                                     ncores = 1,
                                     ...) {

  ##################################################################################################
  ## Checks

  ## Expression Matrix
  if (!is.matrix(object)) stop("'object' must be a matrix.")
  y <- limma::getEAWP(object)
  if (!nrow(y$exprs)) stop("expression matrix has zero rows")

  #	Check weights
  if(!is.null(weights)) {
    stop("weights are not allowed with the phylogenetic regression.")
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

  ##################################################################################################

  tree_model <- get_consensus_tree(y_data, design, phy, model, measurement_error, weights, trim, REML, ncores = ncores, ...)

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
get_consensus_tree <- function(y_data, design, phy, model, measurement_error, weights, trim, REML, ncores, ...) {
  if(!is.null(weights)) stop("weights are not allowed with the phylogenetic regression.")

  if (model == "BM" && !measurement_error) # no parameter to estimate
    return(list(tree = phy,
                params = list(model = "BM",
                              measurement_error = FALSE),
                ddf = rep(nrow(design) - ncol(design), nrow(y_data))))

  all_fits <- fit_all_phylolm(y_data, design, phy, model, measurement_error, weights, REML, ncores, ...)

  get_consensus_tree_model <- switch(model,
                                     BM = get_consensus_tree_BM,
                                     lambda = get_consensus_tree_lambda,
                                     OUfixedRoot = get_consensus_tree_OUfixedRoot)
  params <- get_consensus_tree_model(phy, all_fits, measurement_error, trim)
  params$ddf <- sapply(all_fits, ddf_samples, phylo = phy)

  return(params)
}

#' @title Fit phylolm on all genes
#'
#' @description
#' Fit \code{\link[phylolm]{phylolm}} on all genes.
#'
#' @inheritParams get_C_tree
#'
#' @return A list with all phylolm fits.
#'
#' @keywords internal
#'
fit_all_phylolm <- function(y_data, design, phy, model, measurement_error, weights, REML, ncores, ...) {

  if(!is.null(weights)) stop("weights are not allowed with the phylogenetic regression.")

  alpha_bounds <- getBoundsSelectionStrength(phy)
  min_error <- getMinError(phy)
  lower_bounds <- get_lower_bounds(alpha_bounds, min_error, ...)
  upper_bounds <- get_upper_bounds(alpha_bounds, min_error, ...)
  starting_values <- get_starting_values(alpha_bounds, ...)
  dots_args <- get_dots_args(...)

  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores) # , outfile = ""
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    `%myinfix%` <- `%dopar%`
  } else {
    `%myinfix%` <- `%do%`
  }
  reqpckg <- c("phylolm")

  i <- NULL
  all_fits <- foreach::foreach(i = 1:nrow(y_data), .packages = reqpckg) %myinfix% {

    data_phylolm <- as.data.frame(cbind(y_data[i, ], design))
    colnames(data_phylolm)[1] <- "expr"

    nafun <- function(e) {
      if (grepl("distinguishable", e)) stop(e)
      return(list(sigma2 = NA,
                  sigma2_error = NA,
                  optpar = NA,
                  model = model,
                  n = nrow(data_phylolm),
                  d = ncol(data_phylolm) - 1))
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
    res <- do.call(tmp_fun, dots_args)
    res <- light_phylolm(res)
    rm(data_phylolm)
    return(res)

  }

  return(all_fits)
}

light_phylolm <- function(res) {
  return(list(sigma2 = res$sigma2,
              sigma2_error = res$sigma2_error,
              optpar = res$optpar,
              model = res$model,
              n = res$n,
              d = res$d))
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
  lambda_mean <- tanh(mean_trim(all_lambdas_transform, trim = trim, na.rm = TRUE))
  tree_model <- phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(tree = tree_model,
              params = list(model = "lambda",
                            measurement_error = measurement_error,
                            lambda = lambda_mean,
                            lambda_error = lambda_mean,
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
  all_lambda_transform <- atanh(all_lambda_error)
  # all_lambda_transform <- pracma::logit(all_lambda_error)

  # all_lambdas_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_mean <- tanh(mean_trim(all_lambda_transform, trim = trim, na.rm = TRUE))
  # lambda_mean <- pracma::sigmoid(mean_trim(all_lambda_transform[!is_min_lambda], trim = trim, na.rm = TRUE))

  tree_model <- phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(
    tree = tree_model,
    params = list(model = "BM",
                  measurement_error = measurement_error,
                  lambda_error = lambda_mean,
                  atanh_lambda_error = all_lambda_transform))
  )
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
get_consensus_tree_OUfixedRoot <- function(phy, all_phyfit, measurement_error, trim) {

  # trans_alpha <- function(x) return(atanh(exp(-x)))
  # trans_inv_alpha <- function(x) return(-log(tanh(x)))
  # trans_alpha <- function(x) return(log(x))
  # trans_inv_alpha <- function(x) return(exp(x))
  alpha_bounds <- getBoundsSelectionStrength(phy)
  t_original_tree <- tree_height(phy)
  trans_alpha <- function(alp) {
    atanh(pmax(-1, 1 - rho_prime(alp, t_original_tree)))
  }
  trans_inv_alpha <- function(tt) {
    rho_prime_inv(1 - tanh(tt), t_original_tree, alpha_bounds)
  }

  all_alphas <- sapply(all_phyfit, function(x) x$optpar)
  alpha_na <- !is.finite(all_alphas)
  if (any(alpha_na)) {
    message(paste0("The estimation of the selection strength alpha failed for the following genes, ",
                   "returning a non finite value:\n",
                   paste(which(alpha_na), collapse = ", "),
                   "\nThese genes are excluded from the pool of genes used to build the consensus tree."))
  }
  all_alphas_transform <- trans_alpha(all_alphas)

  if (!measurement_error) {

    alpha_mean <- trans_inv_alpha(mean_trim(all_alphas_transform[!alpha_na], trim = trim, na.rm = TRUE))
    return(list(
      tree = rescale_tree(phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_mean))$tree),
      params = list(model = "OUfixedRoot",
                    measurement_error = measurement_error,
                    alpha = alpha_mean,
                    log_alpha = log(all_alphas),
                    trans_alpha = all_alphas_transform,
                    trans_alpha_fun = "atanh(1-rho)"))
    )
  }

  get_lambda_error_OU <- function(phyfit) {
    if (is.na(phyfit$optpar)) return(NA)
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

  alpha_mean <- trans_inv_alpha(mean_trim(all_alphas_transform, trim = trim, na.rm = TRUE))
  lambda_error_mean <- tanh(mean_trim(all_lambda_error_transform, trim = trim, na.rm = TRUE))

  ## transform tree
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_mean))$tree
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree
  tree_model <- rescale_tree(tree_model)

  return(list(
    tree = tree_model,
    params = list(model = "OUfixedRoot",
                  measurement_error = measurement_error,
                  alpha = alpha_mean,
                  alpha_min = alpha_bounds[1],
                  alpha_max = alpha_bounds[2],
                  trans_alpha = all_alphas_transform,
                  log_alpha = log(all_alphas),
                  trans_alpha_fun = "atanh(1-rho)",
                  lambda_error = lambda_error_mean,
                  atanh_lambda_error = all_lambda_error_transform))
  )
}

#' @title Compute 1 - rho
#'
#' @description
#' rhoprime = 1 - rho = (1 - exp(-2*alpha*t_H)) / (2 * alpha * t_H) is the variance of the OU over the variance of the BM.
#' Taken from Cornuault, 2023, Syst. Biol.
#' It is the fraction of the variance that can be explained by "neutral" BM evolution.
#' When alpha goes to 0, rhoprime goes to 1: trait can be explained by "neutral" BM process.
#' When alpha goes to Inf, rhoprime goes to 0: trait is deterministic.
#'
#' @param alpha the selection strength of the process
#' @param t_tree the total height of the tree
#' @param tol tolerence value for calling alpha = 0.
#'
#' @return Value of rhoprime
#'
#' @keywords internal
#'
rho_prime <- function(alpha, t_tree, tol = .Machine$double.eps) {
  zeros <- alpha <= tol
  nas <- is.na(alpha)
  res <- rep(1, length(alpha))
  res[!zeros & !nas] <- -expm1(-2 * alpha[!zeros & !nas] * t_tree) / (2 * alpha[!zeros & !nas] * t_tree)
  res[nas] <- NA
  return(res)
}

#' @title Compute alpha from rhoprime
#'
#' @description
#' Inverse the rho_prime function.
#'
#' @param y the rhoprime value
#' @param t_tree the total height of the tree
#' @param alpha_bounds lower and upper bounds on alpha values
#'
#' @return Value of alpha
#'
#' @keywords internal
#'
rho_prime_inv <- function(y, t_tree, alpha_bounds) {
  uniroot((function(x) rho_prime(x, t_tree) - y), interval = alpha_bounds, tol = .Machine$double.eps^0.5)$root
}

# #' @title Get OU transformed tree
# #'
# #' @description
# #' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
# #'
# #' @inheritParams get_C_tree
# #'
# #' @return The transformed tree.
# #'
# #' @keywords internal
# #'
# get_consensus_tree_OUrandomRoot <- function(phy, all_phyfit, measurement_error, trim) {
#
#   all_alphas <- sapply(all_phyfit, function(x) x$optpar)
#   all_alphas_transform <- log(all_alphas)
#   # alpha_mean <- exp(mean_trim(all_alphas_transform, trim = trim, na.rm = TRUE))
#
#   if (!measurement_error) {
#
#     alpha_mean <- exp(mean_trim(all_alphas_transform, trim = trim, na.rm = TRUE))
#     return(list(tree = phylolm::transf.branch.lengths(phy, "OUrandomRoot", parameters = list(alpha = alpha_mean))$tree,
#                 params = list(model = "OUrandomRoot",
#                               measurement_error = measurement_error,
#                               alpha = alpha_mean,
#                               log_alpha = all_alphas_transform)))
#   }
#
#   get_lambda_error_OU <- function(phyfit) {
#     tree_model <- phylolm::transf.branch.lengths(phy, "OUrandomRoot",
#                                                  parameters = list(alpha = phyfit$optpar))$tree
#     tilde_t <- tree_height(tree_model) / (2 * phyfit$optpar)
#     lambda_ou_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
#     return(lambda_ou_error)
#   }
#
#   ## consensus lambda error
#   all_lambda_error <- sapply(all_phyfit, get_lambda_error_OU)
#   all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
#   lambda_error_mean <- tanh(mean_trim(all_lambda_error_transform, trim = trim, na.rm = TRUE))
#
#   alpha_mean <- exp(mean_trim(all_alphas_transform, trim = trim, na.rm = TRUE))
#   lambda_error_mean <- tanh(mean_trim(all_lambda_error_transform, trim = trim, na.rm = TRUE))
#
#   ## transform tree
#   tree_model <- phylolm::transf.branch.lengths(phy, "OUrandomRoot", parameters = list(alpha = alpha_mean))$tree
#   tree_model$root.edge <- 0
#   tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree
#   tree_model <- rescale_tree(tree_model)
#
#   return(list(tree = tree_model,
#               params = list(model = "OUrandomRoot",
#                             measurement_error = measurement_error,
#                             alpha = alpha_mean,
#                             lambda_error = lambda_error_mean,
#                             log_alpha = all_alphas_transform,
#                             atanh_lambda_error = all_lambda_error_transform)))
# }

# #' @title Get delta transformed tree
# #'
# #' @description
# #' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
# #'
# #' @inheritParams get_C_tree
# #'
# #' @return The transformed tree.
# #'
# #' @keywords internal
# #'
# get_consensus_tree_delta <- function(phy, all_phyfit, measurement_error, trim) {
#
#   all_deltas <- sapply(all_phyfit, function(x) x$optpar)
#   all_deltas_transform <- log(all_deltas)
#   delta_mean <- exp(mean_trim(all_deltas_transform, trim = trim, na.rm = TRUE))
#
#   if (!measurement_error) {
#
#     return(list(tree = phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = delta_mean))$tree,
#                 params = list(model = "delta",
#                               measurement_error = measurement_error,
#                               delta = delta_mean,
#                               log_delta = all_deltas_transform)))
#   }
#
#   get_lambda_error_delta <- function(phyfit) {
#     tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = phyfit$optpar))$tree
#     tilde_t <- tree_height(tree_model)
#     lambda_delta_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
#     return(lambda_delta_error)
#   }
#
#   all_lambda_error <- sapply(all_phyfit, get_lambda_error_delta)
#   all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
#   lambda_error_mean <- tanh(mean_trim(all_lambda_error_transform, trim = trim, na.rm = TRUE))
#
#   ## transform tree
#   tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = delta_mean))$tree
#   tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree
#   tree_model <- rescale_tree(tree_model)
#
#   return(list(tree = tree_model,
#               params = list(model = "delta",
#                             measurement_error = measurement_error,
#                             delta = delta_mean,
#                             lambda_error = lambda_error_mean,
#                             log_delta = all_deltas_transform,
#                             atanh_lambda_error = all_lambda_error_transform)))
# }

check.consensus_tree <- function(consensus_tree, model, measurement_error) {
  if (consensus_tree$params$model != model) {
    stop(paste0("The consensus tree was computed with the ", consensus_tree$params$model, " model, but the model used is ", model, ". Please check which model is the correct one."))
  }
  if (consensus_tree$params$measurement_error != measurement_error) {
    stop(paste0("The consensus tree was computed with `measurement_error=", consensus_tree$params$measurement_error, "`, but you specified `measurement_error=", measurement_error, "`. Please check which model is the correct one."))
  }
}
