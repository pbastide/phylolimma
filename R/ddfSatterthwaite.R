#' @title Get ddf function
#'
#' @description
#' Get the ddf function
#'
#' @param ddf_method one of "Satterthwaite", "Species", "Samples"
#'
#' @return A function of a phylolm fit and a phylogetic tree
#' returning the ddf.
#'
#' @keywords internal
#'
get_ddf <- function(ddf_method) {
  if (ddf_method == "Satterthwaite") return(ddf_satterthwaite)
  if (ddf_method == "Samples") return(ddf_samples)
  if (ddf_method == "Species") return(ddf_species)
}

#' @title Function for vanilla ddf
#'
#' @param fitlm a phylolm fit
#' @param phylo the corresponding phylogenetic tree
#'
#' @return nsamples - nvariables
#'
#' @keywords internal
#'
ddf_samples <- function(fitlm, phylo) {
  return(fitlm$n - fitlm$d)
}

#' @title Function for species ddf
#'
#' @param fitlm a phylolm fit
#' @param phylo the corresponding phylogenetic tree
#'
#' @return nspecies - nvariables
#'
#' @keywords internal
#'
ddf_species <- function(fitlm, phylo) {
  nspecies <- getSpeciesNumber(phylo)
  return(nspecies - fitlm$d)
}

#' @title Function for Satterthwaite ddf
#'
#' @param fitlm a phylolm fit
#' @param phylo the corresponding phylogenetic tree
#'
#' @return nspecies - nvariables
#'
#' @keywords internal
#'
ddf_satterthwaite <- function(fitlm, phylo) {
  if (!is.null(fitlm$model) && !is.na(fitlm$model) && fitlm$model == "BM" &&
      !is.null(fitlm$sigma2_error) && !is.na(fitlm$sigma2_error) && fitlm$sigma2_error != 0) return(ddf_satterthwaite_BM_error(fitlm, phylo)$ddf[1])
  return(ddf_species(fitlm, phylo))
}

#' @title Satterthwaite for BM with error
#'
#' @description
#' Compute ddf for a \code{\link[phylolm]{phylolm}} fit with the BM and measurement error.
#'
#' @param fit_phylolm a \code{\link[phylolm]{phylolm}} fitted object.
#' Must be fitted with \code{model="BM"} and \code{measurement_error=TRUE}.
#' @param phylo the phylogenetic tre used for the phylolm fit.
#'
#' @return A list, with the approximated degrees of freedom,
#' and the vcov matrix of the parameters.
#'
#' @keywords internal
#'
ddf_satterthwaite_BM_error <- function(fit_phylolm, phylo) {

  X <- fit_phylolm$X
  y <- as.matrix(fit_phylolm$y)
  rownames(y) <- names(fit_phylolm$y)
  yhat <- as.matrix(fit_phylolm$fitted.values)
  rownames(yhat) <- names(fit_phylolm$fitted.values)
  n <- length(phylo$tip.label)
  d <- ncol(X)
  REML <- fit_phylolm$REML

  ## Likelihood
  minusLogLik <- function(pars, y, yhat, X, phy, model) {
    n <- nrow(X)
    d <- ncol(X)
    parameters <- list(sigma2 = exp(pars[1]), sigma2_error = exp(pars[2] - pars[1]))
    phytrans <- phylolm::transf.branch.lengths(phy, model, parameters = parameters)$tree
    comp <- phylolm::three.point.compute(phytrans, P = y - yhat, Q = X)

    if (!REML) {
      n2llh <- as.numeric( n * log(2 * pi) + n * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2) # -2 log-likelihood
    } else {
      # log|X'V^{-1}X|
      ldXX <- determinant(comp$QQ, logarithm = TRUE)$modulus
      n2llh <- as.numeric( (n - d) * log(2 * pi) + (n - d) * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2 + ldXX) # -2 log-likelihood
    }

    return(n2llh / 2)
  }

  optpars <- c(log(fit_phylolm$sigma2), log(fit_phylolm$sigma2_error))

  ## Hessian
  fun <- function(x) {
    return(minusLogLik(x, y, yhat, X, phylo, "BM"))
  }
  J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$sigma2_error))
  A <- compute_inv_hessian(optpars = optpars, fun = fun, grad_trans = J)

  ## Satterthwaite
  K <- ape::vcv(phylo)
  Kd <- diag(diag(K))
  V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
  Vinv <- chol2inv(chol(V))
  ell <- c(0, 1)

  C <- fit_phylolm$vcov
  if (!REML) C <- C * (n - d) / n

  facmat <- C %*% t(X) %*% Vinv
  derfsigma2 <- t(ell) %*% facmat %*% K %*% t(facmat) %*% ell
  derfsigma2_error <- t(ell) %*% facmat %*% Kd %*% t(facmat) %*% ell
  derf <- c(derfsigma2, derfsigma2_error)
  varestim <- t(derf) %*% A %*% derf

  ddf <- 2 * (t(ell) %*% C %*% ell)^2 / varestim

  return(list(ddf = ddf, vcov = A))
}
