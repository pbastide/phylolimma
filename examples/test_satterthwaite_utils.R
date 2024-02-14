bind_star_trees <- function(trees_rep) {
  bind_text <- "("
  tree_text <- sub(";", "", write.tree(trees_rep[1]))
  bind_text <- paste0(bind_text, tree_text)
  for (tt in trees_rep[-1]) {
    tree_text <- sub(";", "", write.tree(tt))
    bind_text <- paste0(bind_text, ",", tree_text)
  }
  bind_text <- paste0(bind_text, ");")
  return(read.tree(text = bind_text))
}

ddf_satterthwaite_sum <- function(fit_phylolm, phylo, REML = FALSE) {

  X <- fit_phylolm$X
  y <- as.matrix(fit_phylolm$y)
  rownames(y) <- names(fit_phylolm$y)
  yhat <- as.matrix(fit_phylolm$fitted.values)
  rownames(yhat) <- names(fit_phylolm$fitted.values)
  n <- length(phylo$tip.label)
  d <- ncol(X)

  # Using the log scale so that parameters are on the entire real line
  optpars <- c(fit_phylolm$sigma2, fit_phylolm$sigma2_error)

  ## Hessian: numerical computation
  K <- vcv(phylo)
  Kd <- diag(diag(K))
  A <- hessianMinusLogLik(optpars, y, yhat, X, K, Kd, REML)
  A <- chol2inv(chol(A))

  ## Gradient of f
  V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
  Vinv <- chol2inv(chol(V))
  ell <- c(0, 1)

  C <- fit_phylolm$vcov
  if (!REML) C <- C * (n - d) / n
  # Cbis <- solve(t(X) %*% Vinv %*% X)
  # all.equal(C, Cbis)

  facmat <- C %*% t(X) %*% Vinv
  derfsigma2 <- t(ell) %*% facmat %*% K %*% t(facmat) %*% ell
  derfsigma2_error <- t(ell) %*% facmat %*% Kd %*% t(facmat) %*% ell
  derf <- c(derfsigma2, derfsigma2_error)

  ## Variance estimation
  varestim <- t(derf) %*% A %*% derf

  ## Satterthwaite
  ddf <- 2 * (t(ell) %*% C %*% ell)^2 / varestim

  return(list(ddf = ddf, vcov = A))
}

hessianMinusLogLik <- function(pars, y, yhat, X, K, Kd, REML) {
  n <- nrow(X)
  d <- ncol(X)
  V <- pars[1] * K + pars[2] * Kd
  Vinv <- chol2inv(chol(V))
  VinvK <- Vinv %*% K
  VinvKsq <- VinvK %*% VinvK
  VinvKd <- Vinv %*% Kd
  VinvKdsq <- VinvKd %*% VinvKd
  VinvKVinvKd <- VinvK %*% VinvKd
  VinvKVinvKdVinv <- VinvKVinvKd %*% Vinv

  dspsp <- - sum(diag(VinvKsq)) + 2 * t(y - yhat) %*% VinvKsq %*% Vinv %*% (y - yhat)
  dsese <- - sum(diag(VinvKdsq)) + 2 * t(y - yhat) %*% VinvKdsq %*% Vinv %*% (y - yhat)
  dspse <- - sum(diag(VinvKVinvKd)) + t(y - yhat) %*% (VinvKVinvKdVinv + t(VinvKVinvKdVinv)) %*% (y - yhat)

  if (REML) {
    tXVinv <- t(X) %*% Vinv
    tXVinvXinv <- solve(tXVinv %*% X)
    tXVinvKVinvX <- tXVinv %*% K %*% t(tXVinv)
    tXVinvKdVinvX <- tXVinv %*% Kd %*% t(tXVinv)

    dspsp <- dspsp - sum(diag(tXVinvXinv %*% tXVinvKVinvX %*% tXVinvXinv %*% tXVinvKVinvX - 2 * tXVinvXinv %*% tXVinv %*% K %*% VinvK %*% t(tXVinv)))
    dsese <- dsese - sum(diag(tXVinvXinv %*% tXVinvKdVinvX %*% tXVinvXinv %*% tXVinvKdVinvX - 2 * tXVinvXinv %*% tXVinv %*% Kd %*% VinvKd %*% t(tXVinv)))
    dspse <- dspse - sum(diag(tXVinvXinv %*% tXVinvKdVinvX %*% tXVinvXinv %*% tXVinvKVinvX - tXVinvXinv %*% tXVinv %*% Kd %*% VinvK %*% t(tXVinv) - tXVinvXinv %*% tXVinv %*% K %*% VinvKd %*% t(tXVinv)))
  }

  return(matrix(c(dspsp, dspse, dspse, dsese), 2, 2) / 2)
}

ddf_satterthwaite_sum_approx <- function(fit_phylolm, phylo, REML = FALSE) {

  X <- fit_phylolm$X
  y <- as.matrix(fit_phylolm$y)
  rownames(y) <- names(fit_phylolm$y)
  yhat <- as.matrix(fit_phylolm$fitted.values)
  rownames(yhat) <- names(fit_phylolm$fitted.values)
  n <- length(phylo$tip.label)
  d <- ncol(X)

  ## Likelihood function
  minusLogLik <- function(pars, y, yhat, X, phy, model) {
    n <- nrow(X)
    d <- ncol(X)
    parameters <- list(sigma2 = exp(pars[1]), sigma2_error = exp(pars[2] - pars[1]))
    phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
    comp <- three.point.compute(phytrans, P = y - yhat, Q = X)

    if (!REML) {
      n2llh <- as.numeric( n * log(2 * pi) + n * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2) # -2 log-likelihood
    } else {
      # log|X'V^{-1}X|
      ldXX <- determinant(comp$QQ, logarithm = TRUE)$modulus
      n2llh <- as.numeric( (n - d) * log(2 * pi) + (n - d) * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2 + ldXX) # -2 log-likelihood
    }

    return(n2llh / 2)
  }

  # Using the log scale so that parameters are on the entire real line
  optpars <- c(log(fit_phylolm$sigma2), log(fit_phylolm$sigma2_error))

  # all.equal(minusLogLik(optpars, y, yhat, X, phylo, "BM"),
  #           -fit_phylolm$logLik)

  ## Hessian: numerical computation
  fun <- function(x) {
    return(minusLogLik(x, y, yhat, X, phylo, "BM"))
  }
  J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$sigma2_error))
  A <- compute_hessian(optpars = optpars, fun = fun, grad_trans = J)

  ## Gradient of f
  K <- vcv(phylo)
  Kd <- diag(diag(K))
  V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
  Vinv <- chol2inv(chol(V))
  ell <- c(0, 1)

  C <- fit_phylolm$vcov
  if (!REML) C <- C * (n - d) / n
  # Cbis <- solve(t(X) %*% Vinv %*% X)
  # all.equal(C, Cbis)

  facmat <- C %*% t(X) %*% Vinv
  derfsigma2 <- t(ell) %*% facmat %*% K %*% t(facmat) %*% ell
  derfsigma2_error <- t(ell) %*% facmat %*% Kd %*% t(facmat) %*% ell
  derf <- c(derfsigma2, derfsigma2_error)

  ## Variance estimation
  varestim <- t(derf) %*% A %*% derf

  ## Satterthwaite
  ddf <- 2 * (t(ell) %*% C %*% ell)^2 / varestim

  return(list(ddf = ddf, vcov = A))
}

# Adapted from lmerTest
# https://github.com/runehaubo/lmerTestR/blob/35dc5885205d709cdc395b369b08ca2b7273cb78/R/lmer.R#L173
compute_hessian <- function(optpars, fun, grad_trans, tol = 1e-8, ...) {
  # Compute Hessian:
  h <- numDeriv::hessian(func = fun, x = optpars, ...)
  # back transformation of parameters
  h <- t(grad_trans) %*% h %*% grad_trans
  # Eigen decompose the Hessian:
  eig_h <- eigen(h, symmetric=TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  zero <- evals > -tol & evals < tol
  if(sum(neg) > 0) { # negative eigenvalues
    eval_chr <- if(sum(neg) > 1) "eigenvalues" else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf("Model failed to converge with %d negative %s: %s",
                    sum(neg), eval_chr, evals_num), call.=FALSE)
  }
  # Note: we warn about negative AND zero eigenvalues:
  if(sum(zero) > 0) { # some eigenvalues are zero
    eval_chr <- if(sum(zero) > 1) "eigenvalues" else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf("Model may not have converged with %d %s close to zero: %s",
                    sum(zero), eval_chr, evals_num))
  }
  # Compute vcov(varpar):
  pos <- eig_h$values > tol
  q <- sum(pos)
  # Using the Moore-Penrose generalized inverse for h:
  h_inv <- with(eig_h, {
    vectors[, pos, drop=FALSE] %*% diag(1/values[pos], nrow=q) %*%
      t(vectors[, pos, drop=FALSE]) })
  return(h_inv)
}
