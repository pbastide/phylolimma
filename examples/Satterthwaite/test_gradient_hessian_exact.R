minusLogLik <- function(pars, y, yhat, X, phy, model, REML) {
  n <- nrow(X)
  d <- ncol(X)
  parameters <- list(sigma2 = pars[1], sigma2_error = pars[2] / pars[1])
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

gradientMinusLogLik <- function(pars, y, yhat, X, K, Kd, REML) {
  n <- nrow(X)
  d <- ncol(X)
  V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
  Vinv <- chol2inv(chol(V))


  dsigmaphylo <- sum(diag(Vinv %*% K)) - t(y - yhat) %*% Vinv %*% K %*% Vinv %*% (y - yhat)
  dsigmaerror <- sum(diag(Vinv %*% Kd)) - t(y - yhat) %*% Vinv %*% Kd %*% Vinv %*% (y - yhat)

  if (REML) {
    tXVinv <- t(X) %*% Vinv
    tXVinvXinv <- solve(tXVinv %*% X)
    dsigmaphylo <- dsigmaphylo - sum(diag(tXVinvXinv %*% tXVinv %*% K %*% t(tXVinv)))
    dsigmaerror <- dsigmaerror - sum(diag(tXVinvXinv %*% tXVinv %*% Kd %*% t(tXVinv)))
  }

  return(c(dsigmaphylo, dsigmaerror) / 2)
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

################################################################################
## Test
################################################################################

library(lmerTest)
library(phylolm)

set.seed(12891289)

## Star Tree with replicates
ntips <- 4
reps = c(3, 4, 2, 5) # replicates per species
ntot = sum(reps)
eps <- 1e-16 # almost zero
trees_rep <- vector(mode = "list", length = ntips)
for (r in 1:ntips) {
  trees_rep[[r]] <- stree(reps[r]) # star tree for replicates
  trees_rep[[r]]$edge.length <- rep(eps, nrow(trees_rep[[r]]$edge)) # almost zero distance
  trees_rep[[r]]$root.edge <- 1 # root at length 1 from tips
  trees_rep[[r]]$tip.label <- paste0("t", r, "_", 1:reps[r])
}
tree_rep <- bind_star_trees(trees_rep)
plot(tree_rep)

## Group Design
traits <- data.frame(species = sub("_[0-9]", "", tree_rep$tip.label),
                     id = tree_rep$tip.label)
rownames(traits) <- traits$id
design <- paste0("t", c(1, 2))
traits$design <- traits$species %in% design

plot(tree_rep)
tiplabels(pch = 21, col = as.factor(traits$design + 0), bg = as.factor(traits$design + 0))

## Traits
sigma2_phylo <- 0.8
sigma2_error <- 0.4
traits$g1 <- phylolm::rTrait(n = 1, tree_rep, model = "BM",
                             parameters = list(ancestral.state = 0,
                                               sigma2 = sigma2_phylo))
traits$g1 <- traits$g1 + rnorm(length(traits$g1), mean = 0, sd = sqrt(sigma2_error))
# Null model : no difference between groups
traits$g1[traits$design] <- traits$g1[traits$design] + 0

## phylolm
REML <- TRUE
fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep,
                                model = "BM", measurement_error = TRUE,
                                REML = REML)

X <- fit_phylolm$X
y <- as.matrix(fit_phylolm$y)
rownames(y) <- names(fit_phylolm$y)
yhat <- as.matrix(fit_phylolm$fitted.values)
rownames(yhat) <- names(fit_phylolm$fitted.values)
n <- length(tree_rep$tip.label)
d <- ncol(X)
K <- vcv(tree_rep)
Kd <- diag(diag(K))

optpars <- c(fit_phylolm$sigma2, fit_phylolm$sigma2_error)
fun <- function(x) {
  return(minusLogLik(x, y, yhat, X, tree_rep, "BM", REML))
}

numDeriv::grad(func = fun, x = optpars)
gradientMinusLogLik(optpars, y, yhat, X, K, Kd, REML = REML)

numDeriv::hessian(func = fun, x = optpars)
hessianMinusLogLik(optpars, y, yhat, X, K, Kd, REML = REML)
