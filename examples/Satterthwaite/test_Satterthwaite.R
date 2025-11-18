library(phylolm)
## Tree
set.seed(1289)
ntips <- 10
tree <- ape::rphylo(ntips, 0.1, 0)
tree$edge.length <- tree$edge.length / vcv(tree)[1, 1]

plot(tree)
nodelabels()

r <- 3
traits <- data.frame(species = rep(tree$tip.label, r))
traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
rownames(traits) <- traits$id

tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)

traits <- traits[match(tree_rep$tip.label, traits$id), ]

# design <- paste0("t", c(4, 11, 12, 13, 15, 7, 14, 6, 10, 18, 13, 16))
design <- paste0("t", c(3, 8))
# design <- paste0("t", c(3, 8, 5, 9, 4))
# design <- paste0("t", c(3, 5, 4, 1, 10))
# design <- paste0("t", c(4, 11))
# design <- tree$tip.label[phytools::getDescendants(tree, node = 107)]
design <- design[!is.na(design)]
design <- traits$species %in% design
traits$design <- as.factor(design + 0)

plot(tree_rep)
tiplabels(pch = 21, col = traits$design, bg = traits$design)

##############################################################################
##############################################################################
### Satterthwaite sum
##############################################################################
##############################################################################
Nrep <- 100
pvanilla <- psatt <- psatt_lambda <- pbis <- dsatt <- dsatt_lambda <- NULL
effect <- 0
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
  traits$g1 <- sim + rnorm(length(sim), 0, sd = sqrt(0.2))
  traits$g1[traits$design] <- traits$g1[traits$design] + effect

  n <- length(tree_rep$tip.label)
  d <- 2

  ## phylolm BM
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE)
  summary(fit_phylolm)

  res_satt <- ddf_satterthwaite_sum(fit_phylolm, tree_rep)
  df <- res_satt$ddf
  dsatt <- c(dsatt, df)

  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
  pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = n-d))
  psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
  pbis <- c(pbis, 2 * pt(-abs(tval), df = ntips - d))

  ## phylolm lambda
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
                                         measurement_error = FALSE,
                                         lower.bound = list(lambda = 1e-16),
                                         upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))))

  res_satt_lambda <- ddf_satterthwaite_lambda(fit_phylolm_lambda, tree_rep)
  df_lambda <- res_satt_lambda$ddf
  dsatt_lambda <- c(dsatt_lambda, df_lambda)

  tval <- summary(fit_phylolm_lambda)$coefficients[2, "t.value"]
  psatt_lambda <- c(psatt_lambda, 2 * pt(-abs(tval), df = df_lambda))
#
#   ## GLS
#   fit_gls <- nlme::gls(g1 ~ design, traits, correlation = corPagel(1, tree_rep, form = ~id))
#
#   res_satt_lambda <- ddf_satterthwaite_lambda(fit_gls, tree_rep)
#   df_lambda <- res_satt_lambda$ddf
#   dsatt_lambda <- c(dsatt_lambda, df_lambda)
#
#   tval <- summary(fit_phylolm_lambda)$coefficients[2, "t.value"]
#   psatt_lambda <- c(psatt_lambda, 2 * pt(-abs(tval), df = df_lambda))

}
hist(pvanilla)
hist(psatt)
hist(pbis)
hist(psatt_lambda)
plot(sort(pvanilla))
points(sort(psatt), col = "red")
points(sort(psatt_lambda), col = "orange")
points(sort(pbis), col = "green")
abline(a = 0, b = 1/Nrep)
plot(dsatt, ylim = c(0, fit_phylolm$n), col = "red")
points(dsatt_lambda, col = "orange")
abline(a = fit_phylolm$n-2, b = 0)
abline(a = ntips-2, b = 0, col = "green")
mean(psatt <= 0.05)
mean(psatt_lambda <= 0.05)
mean(pvanilla <= 0.05)
mean(pbis <= 0.05)

##############################################################################
##############################################################################
### Satterthwaite REML
##############################################################################
##############################################################################
Nrep <- 100
pvanilla <- psatt <- psatt_lambda <- pbis <- dsatt <- dsatt_lambda <- NULL
effect <- 0
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
  traits$g1 <- sim + rnorm(length(sim), 0, sd = sqrt(0.2))
  traits$g1[traits$design] <- traits$g1[traits$design] + effect

  n <- length(tree_rep$tip.label)
  d <- 2

  ## phylolm BM
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE, REML = TRUE)
  summary(fit_phylolm)

  res_satt <- ddf_satterthwaite_sum(fit_phylolm, tree_rep, REML = TRUE)
  df <- res_satt$ddf
  dsatt <- c(dsatt, df)

  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
  pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = n-d))
  psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
  pbis <- c(pbis, 2 * pt(-abs(tval), df = ntips - d))

  ## phylolm lambda
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
                                         measurement_error = FALSE,
                                         lower.bound = list(lambda = 1e-16),
                                         upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))), REML = TRUE)

  res_satt_lambda <- ddf_satterthwaite_lambda(fit_phylolm_lambda, tree_rep, REML = TRUE)
  df_lambda <- res_satt_lambda$ddf
  dsatt_lambda <- c(dsatt_lambda, df_lambda)

  tval <- summary(fit_phylolm_lambda)$coefficients[2, "t.value"]
  psatt_lambda <- c(psatt_lambda, 2 * pt(-abs(tval), df = df_lambda))

}
hist(pvanilla)
hist(psatt)
hist(pbis)
hist(psatt_lambda)
plot(sort(pvanilla))
points(sort(psatt), col = "red")
points(sort(psatt_lambda), col = "orange")
points(sort(pbis), col = "green")
abline(a = 0, b = 1/Nrep)
plot(dsatt, ylim = c(0, fit_phylolm$n), col = "red")
points(dsatt_lambda, col = "orange")
abline(a = fit_phylolm$n-2, b = 0)
abline(a = ntips-2, b = 0, col = "green")
mean(psatt <= 0.05)
mean(psatt_lambda <= 0.05)
mean(pvanilla <= 0.05)
mean(pbis <= 0.05)

##############################################################################
##############################################################################
## No phylo
##############################################################################
##############################################################################

Nrep <- 100
plm <- NULL
effect <- 0
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  r <- 2
  n <- length(tree$tip.label) * r
  traits <- data.frame(species = rep(tree$tip.label, r))
  traits$g1 <- rnorm(n, 0, 0.2)

  design <- paste0("t", c(1, 3))
  traits$design <- traits$species %in% design
  traits$g1[traits$design] <- traits$g1[traits$design] + effect
  traits$design <- as.factor(traits$design + 0)

  ### trait

  ## phylolm
  fit_lm <- lm(g1 ~ design, traits)
  summary(fit_lm)

  tval <- summary(fit_lm)$coefficients[2, "t value"]
  plm <- c(plm, 2 * pt(-abs(tval), df = n - 2))
}
hist(plm)
plot(sort(plm))
abline(a = 0, b = 1/Nrep)

##############################################################################
##############################################################################
## Simple BM
##############################################################################
##############################################################################

Nrep <- 100
pBM <- NULL
effect <- 0
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  r <- 1
  traits <- data.frame(species = rep(tree$tip.label, r))
  traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
  rownames(traits) <- traits$id

  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)

  traits <- traits[match(tree_rep$tip.label, traits$id), ]

  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 1))
  traits$g1 <- sim + rnorm(length(sim), 0, 0)

  # design <- c("t4", "t11", "t12", "t13", "t15", "t7", "t14", "t6", "t10", "t18", "t13", "t16")
  design <- paste0("t", c(3, 8))
  traits$design <- traits$species %in% design
  traits$g1[traits$design] <- traits$g1[traits$design] + effect
  traits$design <- as.factor(traits$design + 0)

  ## Plot
  plot(tree_rep)
  tiplabels(pch = 21, col = traits$design, bg = traits$design)

  ##############################################################################
  ### BM

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = FALSE)
  summary(fit_phylolm)

  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
  pBM <- c(pBM, 2 * pt(-abs(tval), df = n-d))
}
hist(pBM)
plot(sort(pBM))
abline(a = 0, b = 1/Nrep)

# ##############################################################################
# ##############################################################################
# ### Satterthwaite BM error
# ##############################################################################
# ##############################################################################
# Nrep <- 100
# pvanilla <- psatt <- dsatt <- NULL
# effect <- 0
# set.seed(1289)
# for (rep in 1:Nrep) {
#
#   ## Simulate data
#   sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
#   traits$g1 <- sim + rnorm(length(sim), 0, sd = sqrt(0.2))
#   traits$g1[traits$design] <- traits$g1[traits$design] + effect
#
#   ### BM
#
#   ## phylolm
#   fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE)
#   summary(fit_phylolm)
#
#   ## Likelihood
#   X <- fit_phylolm$X
#   y <- as.matrix(traits$g1)
#   rownames(y) <- traits$id
#
#   minusLogLik <- function(pars, y, X, phy, model) {
#     parameters <- list(sigma2_error = exp(pars[1]))
#     phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
#     n <- length(phytrans$tip.label)
#
#     comp <- three.point.compute(phytrans, P = y, Q = X)
#
#     invXX <- solve(comp$QQ)
#     betahat <- invXX %*% comp$QP
#     sigma2hat <- as.numeric((comp$PP - 2 * t(betahat) %*% comp$QP + t(betahat) %*% comp$QQ %*% betahat) / n)
#     if (sigma2hat<0) {
#       resdl <- X %*% betahat - y
#       compyy <- three.point.compute(phytrans, P = resdl, Q = X)
#       sigma2hat <- compyy$PP / n
#     }
#     n2llh <- as.numeric( n * log(2 * pi) + n + n * log(sigma2hat) + comp$logd) # -2 log-likelihood
#     return(n2llh / 2)
#   }
#
#   optpars <- c(log(fit_phylolm$sigma2_error) - log(fit_phylolm$sigma2))
#
#   all.equal(minusLogLik(optpars, y, X, tree_rep, "BM"),
#             -fit_phylolm$logLik)
#
#   ## Hessian
#   fun <- function(x) {
#     return(minusLogLik(x, y, X, tree_rep, "lambda"))
#   }
#   approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/5))
#   A <- 1 / approxHessian$Hessian[1, 1] / (fit_phylolm$sigma2 / fit_phylolm$sigma2_error)^2
#
#   # approxHessian <- pracma::hessian(f = minusLogLik, x0 = optpars,
#   #                                  y = y, X = X, phy = tree_rep, model = "BM")
#   # A <- 1 / approxHessian[1, 1] / (fit_phylolm$sigma2 / fit_phylolm$sigma2_error)^2
#
#   ## Satterthwaite
#   n <- length(tree_rep$tip.label)
#   d <- 2
#   K <- vcv(tree_rep)
#   I <- diag(rep(1, n))
#   V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * I
#   Vinv <- solve(V)
#   gamma <- fit_phylolm$sigma2_error / fit_phylolm$sigma2
#   W <- K + gamma * I
#   Winv <- solve(W)
#   ell <- c(0, 1)
#
#   C <- fit_phylolm$vcov * (n - 2) / n
#   Cbis <- solve(t(X) %*% Vinv %*% X)
#   all.equal(C, Cbis)
#
#   D <- solve(t(X) %*% Winv %*% X)
#   all.equal(D, C / fit_phylolm$sigma2)
#
#   facmat <- D %*% t(X) %*% Winv
#   derCgamma <- facmat %*% I %*% t(facmat)
#   derfgamma <- t(ell) %*% derCgamma %*% ell
#
#   dfsigerrinv <- derfgamma^2 * A / 2 / (t(ell) %*% D %*% ell)^2
#
#   df <- (1 / (n - d) + dfsigerrinv)^{-1}
#   # if (is.infinite(A)) df <- 0.1
#   dsatt <- c(dsatt, df)
#
#   tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
#   pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = n-d))
#   psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
# }
# hist(pvanilla)
# hist(psatt)
# plot(sort(pvanilla))
# points(sort(psatt), col = "red")
# abline(a = 0, b = 1/Nrep)
# plot((n-2) - dsatt)
#
# ##############################################################################
# ##############################################################################
# ### Satterthwaite lambda
# ##############################################################################
# ##############################################################################
# Nrep <- 100
# pvanilla <- psatt <- dsatt <- NULL
# effect <- 0
# set.seed(1289)
# for (rep in 1:Nrep) {
#
#   ## Simulate data
#   r <- 3
#   traits <- data.frame(species = rep(tree$tip.label, r))
#   traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
#   rownames(traits) <- traits$id
#
#   tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)
#
#   traits <- traits[match(tree_rep$tip.label, traits$id), ]
#
#   sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
#   traits$g1 <- sim + rnorm(length(sim), 0, sd = sqrt(0.2))
#
#   design <- paste0("t", c(4, 11, 12, 13, 15, 7, 14, 6, 10, 18, 13, 16))
#   # design <- paste0("t", c(1, 3))
#   traits$design <- traits$species %in% design
#   traits$g1[traits$design] <- traits$g1[traits$design] + effect
#   traits$design <- as.factor(traits$design + 0)
#
#   ## Plot
#   plot(tree_rep)
#   tiplabels(pch = 21, col = traits$design, bg = traits$design)
#
#   ### BM
#
#   ## phylolm
#   fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda", lower.bound = list(lambda = 1e-16))
#   summary(fit_phylolm)
#
#   ## Likelihood
#   X <- fit_phylolm$X
#   y <- as.matrix(traits$g1)
#   rownames(y) <- traits$id
#
#   # minusLogLik <- function(pars, y, X, phy, model) {
#   #   parameters <- list(sigma2 = pars[1], lambda = pars[2])
#   #   phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
#   #   n <- length(phytrans$tip.label)
#   #   comp <- three.point.compute(phytrans, P = y, Q = X)
#   #   # if (parameters$sigma2 <= 0) return(1000)
#   #   n2llh <- as.numeric( n * log(2 * pi) + n + n * log(parameters$sigma2) + comp$logd) # -2 log-likelihood
#   #   return(n2llh / 2)
#   # }
#
#   minusLogLik <- function(pars, y, X, phy, model) {
#     parameters <- list(lambda = exp(pars[1]))
#     phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
#     n <- length(phytrans$tip.label)
#
#     comp <- three.point.compute(phytrans, P = y, Q = X)
#
#     invXX <- solve(comp$QQ)
#     betahat <- invXX %*% comp$QP
#     sigma2hat <- as.numeric((comp$PP - 2 * t(betahat) %*% comp$QP + t(betahat) %*% comp$QQ %*% betahat) / n)
#     if (sigma2hat<0) {
#       resdl <- X %*% betahat - y
#       compyy <- three.point.compute(phytrans, P = resdl, Q = X)
#       sigma2hat <- compyy$PP / n
#     }
#     n2llh <- as.numeric( n * log(2 * pi) + n + n * log(sigma2hat) + comp$logd) # -2 log-likelihood
#     return(n2llh / 2)
#   }
#
#   optpars <- c(log(fit_phylolm$optpar))
#
#   all.equal(minusLogLik(optpars, y, X, tree_rep, "lambda"),
#             -fit_phylolm$logLik)
#
#   ## Hessian
#   fun <- function(x) {
#     return(minusLogLik(x, y, X, tree_rep, "lambda"))
#   }
#   approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/4))
#   A <- 1 / approxHessian$Hessian[1, 1] / (1/fit_phylolm$optpar)^2
#
#   # approxHessian <- pracma::hessian(f = minusLogLik, x0 = optpars,
#   #                                  y = y, X = X, phy = tree_rep, model = "lambda", h = .Machine$double.eps^(1/4))
#   # A <- 1 / approxHessian[1, 1] / (1/fit_phylolm$optpar)^2
#
#   ## Satterthwaite
#   n <- length(tree_rep$tip.label)
#   d <- 2
#   K <- vcv(tree_rep)
#   Kd <- diag(diag(K))
#   lambda <- fit_phylolm$optpar
#   W <- lambda * K + (1 - lambda) * Kd
#   Winv <- solve(W)
#   ell <- c(0, 1)
#
#   D <- solve(t(X) %*% Winv %*% X)
#   all.equal(D, fit_phylolm$vcov * (n - 2) / n / fit_phylolm$sigma2)
#
#   facmat <- D %*% t(X) %*% Winv
#   derWgamma <- K - Kd
#   derDgamma <- facmat %*% derWgamma %*% t(facmat)
#   derfgamma <- t(ell) %*% derDgamma %*% ell
#
#   dfsigerrinv <- derfgamma^2 * A / 2 / (t(ell) %*% D %*% ell)^2
#
#   df <- (1 / (n - d) + dfsigerrinv)^{-1}
#   # if (is.infinite(A)) df <- 0.1
#   dsatt <- c(dsatt, df)
#
#   tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
#   pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = n-d))
#   psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
# }
# hist(pvanilla)
# hist(psatt)
# plot(sort(pvanilla))
# points(sort(psatt), col = "red")
# abline(a = 0, b = 1/Nrep)
# plot((n-2) - dsatt)
