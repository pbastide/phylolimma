
ddf_satterthwaite_BM <- function(fit_phylolm, phylo, REML) {

  X <- fit_phylolm$X
  y <- as.matrix(fit_phylolm$y)
  rownames(y) <- names(fit_phylolm$y)

  ## Likelihood
  minusLogLik <- function(pars, y, X, phy, model) {
    parameters <- list(sigma2_error = exp(pars[1]))
    # parameters <- list(sigma2_error = pars[1])
    phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
    n <- nrow(X)
    d <- ncol(X)

    comp <- three.point.compute(phytrans, P = y, Q = X)

    invXX <- solve(comp$QQ)
    betahat <- invXX %*% comp$QP
    sigma2hat <- as.numeric((comp$PP - 2 * t(betahat) %*% comp$QP + t(betahat) %*% comp$QQ %*% betahat) / n)
    if (REML) sigma2hat <- sigma2hat * n / (n - d)
    if (sigma2hat<0) {
      resdl <- X %*% betahat - y
      compyy <- three.point.compute(phytrans, P = resdl, Q = X)
      sigma2hat <- compyy$PP / n
      if (REML) sigma2hat <- sigma2hat * n / (n - d)
    }
    if (!REML) {
      n2llh <- as.numeric( n * log(2 * pi) + n * log(sigma2hat) + n + comp$logd) # -2 log-likelihood
    } else {
      # log|X'V^{-1}X|
      ldXX <- determinant(comp$QQ, logarithm = TRUE)$modulus
      n2llh <- as.numeric( (n - d) * log(2 * pi) + (n - d) * log(sigma2hat) + n - d + comp$logd + ldXX) # -2 log-likelihood
    }
    return(n2llh / 2)
  }

  optpars <- c(log(fit_phylolm$sigma2_error) - log(fit_phylolm$sigma2))
  # optpars <- c(fit_phylolm$sigma2_error / fit_phylolm$sigma2)

  # all.equal(minusLogLik(optpars, y, X, phylo, "BM"),
  #           -fit_phylolm$logLik)

  ## Hessian
  fun <- function(x) {
    return(minusLogLik(x, y, X, phylo, "BM"))
  }
  approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/5))
  A <- 1 / approxHessian$Hessian[1, 1] / (fit_phylolm$sigma2 / fit_phylolm$sigma2_error)^2

  # approxHessian <- pracma::hessian(f = minusLogLik, x0 = fit_phylolm$sigma2_error / fit_phylolm$sigma2,
  #                                  y = y, X = X, phy = phylo, model = "BM")
  # A <- 1 / approxHessian[1, 1]

  ## Satterthwaite
  n <- length(phylo$tip.label)
  d <- 2
  K <- vcv(phylo)
  Kd <- diag(diag(K))
  V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
  Vinv <- chol2inv(chol(V))
  gamma <- fit_phylolm$sigma2_error / fit_phylolm$sigma2
  W <- K + gamma * Kd
  Winv <- chol2inv(chol(W))
  ell <- c(0, 1)

  # C <- fit_phylolm$vcov * (n - 2) / n
  # Cbis <- solve(t(X) %*% Vinv %*% X)
  # all.equal(C, Cbis)

  D <- solve(t(X) %*% Winv %*% X)
  # all.equal(D, C / fit_phylolm$sigma2)

  facmat <- D %*% t(X) %*% Winv
  derDgamma <- facmat %*% Kd %*% t(facmat)
  derfgamma <- t(ell) %*% derDgamma %*% ell

  dfsigerrinv <- derfgamma^2 * A / 2 / (t(ell) %*% D %*% ell)^2 * (1 + 2 / (n - d))

  df <- (1 / (n - d) + dfsigerrinv)^{-1}
  return(list(df = df, vcov = A))
}

ddf_satterthwaite_sum <- function(fit_phylolm, phylo, REML = FALSE) {

  X <- fit_phylolm$X
  y <- as.matrix(fit_phylolm$y)
  rownames(y) <- names(fit_phylolm$y)
  yhat <- as.matrix(fit_phylolm$fitted.values)
  rownames(yhat) <- names(fit_phylolm$fitted.values)
  n <- length(phylo$tip.label)
  d <- ncol(X)

  ## Likelihood
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


  optpars <- c(log(fit_phylolm$sigma2), log(fit_phylolm$sigma2_error))

  # all.equal(minusLogLik(optpars, y, yhat, X, phylo, "BM"),
  #           -fit_phylolm$logLik)

  ## Hessian
  fun <- function(x) {
    return(minusLogLik(x, y, yhat, X, phylo, "BM"))
  }
  J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$sigma2_error))
  A <- compute_hessian(optpars = optpars, fun = fun, grad_trans = J)
  # approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/5))
  # H <- approxHessian$Hessian
  # J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$sigma2_error))
  # H <- t(J) %*% approxHessian$Hessian %*% J
  # A <- solve(H)

  # approxHessian <- pracma::hessian(f = minusLogLik, x0 = fit_phylolm$sigma2_error / fit_phylolm$sigma2,
  #                                  y = y, X = X, phy = phylo, model = "BM")
  # A <- 1 / approxHessian[1, 1]

  ## Satterthwaite
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
  varestim <- t(derf) %*% A %*% derf

  ddf <- 2 * (t(ell) %*% C %*% ell)^2 / varestim

  return(list(ddf = ddf, vcov = A))
}

ddf_satterthwaite_lambda <- function(fit_phylolm, phylo, REML = FALSE) {

  X <- fit_phylolm$X
  y <- as.matrix(fit_phylolm$y)
  rownames(y) <- names(fit_phylolm$y)
  yhat <- as.matrix(fit_phylolm$fitted.values)
  rownames(yhat) <- names(fit_phylolm$fitted.values)
  n <- length(phylo$tip.label)
  d <- ncol(X)

  ## Likelihood
  minusLogLik <- function(pars, y, yhat, X, phy, model) {
    n <- nrow(X)
    d <- ncol(X)
    parameters <- list(sigma2 = exp(pars[1]), lambda = exp(pars[2]))
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

  optpars <- c(log(fit_phylolm$sigma2), log(fit_phylolm$optpar))

  # all.equal(minusLogLik(optpars, y, yhat, X, phylo, "lambda"),
  #           -fit_phylolm$logLik)

  ## Hessian
  fun <- function(x) {
    return(minusLogLik(x, y, yhat, X, phylo, "lambda"))
  }
  J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$optpar))
  A <- compute_hessian(optpars = optpars, fun = fun, grad_trans = J)
  # approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/3))
  # H <- approxHessian$Hessian
  # J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$optpar))
  # H <- t(J) %*% approxHessian$Hessian %*% J
  # A <- solve(H)

  ## Satterthwaite
  K <- vcv(phylo)
  Kd <- diag(diag(K))
  lambda <- fit_phylolm$optpar
  V <- fit_phylolm$sigma2 * (lambda * K + (1 - lambda) * Kd)
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
  varestim <- t(derf) %*% A %*% derf

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
