library(phylolm)

## Tree Chen small
tree_rep <- read.tree(text = "(((opossum.2:1.359619307e-14,opossum.1:0):1.359619307e-14,opossum.0:0):1,(((cow.2:1.359619307e-14,cow.1:0):1.359619307e-14,cow.0:0):0.6539768865,(((rat.2:1.359619307e-14,rat.1:0):1.359619307e-14,rat.0:0):0.6247450714,(human.1:1.359619307e-14,human.0:0):0.6247450714):0.02923181509):0.3460231135);")
tree_rep$edge.length <- tree_rep$edge.length / vcv(tree_rep)[1,1]
species <- tree_rep$tip.label
design <- paste0("rat.", 0:2)
design <- species %in% design
plot(tree_rep)
tiplabels(pch = 21, col = as.factor(design + 0), bg = as.factor(design + 0))
ntips <- 4

# ## Tree Chen full
# tree_rep <- read.tree(text = "(((opossum.2:2e-12,opossum.1:0):2e-12,opossum.0:0):147.1,(armadillo.0:98.7,((((cow.2:2e-12,cow.1:0):2e-12,cow.0:0):85.4,(ferret.0:55.9,dog.0:55.9):29.5):10.8,((rabbit.0:89,(((rat.2:2e-12,rat.1:0):2e-12,rat.0:0):29.4,((((musCaroli.3:2e-12,musCaroli.2:0):2e-12,musCaroli.1:0):2e-12,musCaroli.0:0):4.1,(((((musMusculus.5:2e-12,musMusculus.4:0):2e-12,musMusculus.3:0):2e-12,musMusculus.2:0):2e-12,musMusculus.1:0):2e-12,musMusculus.0:0):4.1,((musSpretus.2:2e-12,musSpretus.1:0):2e-12,musSpretus.0:0):4.1):25.3):59.6):2.9,((marmoset.1:2e-12,marmoset.0:0):52.4,(((((rhesus.4:2e-12,rhesus.3:0):2e-12,rhesus.2:0):2e-12,rhesus.1:0):2e-12,rhesus.0:0):35.3,((orangutan.1:2e-12,orangutan.0:0):18.4,((gorilla.1:2e-12,gorilla.0:0):11.5,((human.1:2e-12,human.0:0):8.8,((bonobo.1:2e-12,bonobo.0:0):4.5,(chimp.1:2e-12,chimp.0:0):4.5):4.3):2.7):6.9):16.9):17.1):39.5):4.3):2.5):48.4);")
# tree_rep$edge.length <- tree_rep$edge.length / vcv(tree_rep)[1,1]
# species <- tree_rep$tip.label
# design <- c(paste0("rat.", 0:8),
#             paste0("musMusculus.", 0:8),
#             paste0("musCaroli.", 0:8),
#             paste0("musSpretus.", 0:8))
# design <- species %in% design
# plot(tree_rep)
# tiplabels(pch = 21, col = as.factor(design + 0), bg = as.factor(design + 0))
# ntips <- 17

##############################################################################
##############################################################################
### Satterthwaite BM
##############################################################################
##############################################################################
Nrep <- 100
pvanilla <- psatt <- psatt_lambda <- pbis <- dsatt <- dsatt_lambda <- NULL
effect <- 0
sigphylo <- 0.8
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  traits <- data.frame(species = species, id = species)

  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = sigphylo))

  traits <- data.frame(species = species,
                       id = species,
                       g1 = sim + rnorm(length(sim), 0, sd = sqrt(1 - sigphylo)),
                       design = as.factor(design + 0))
  rownames(traits) <- species
  traits$g1[design] <- traits$g1[design] + effect

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM",
                                  measurement_error = TRUE,
                                  lower.bound = list(sigma2_error = getMinError(tree_rep)))
  # summary(fit_phylolm)
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
                                         measurement_error = FALSE,
                                         lower.bound = list(lambda = 1e-16),
                                         upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))))
  # fit_gls <- nlme::gls(g1 ~ design, traits, corPagel(0.5, tree_rep, form = ~species))

  res_satt <- ddf_satterthwaite_sum(fit_phylolm, tree_rep)
  df <- res_satt$ddf
  dsatt <- c(dsatt, df)

  res_satt_lambda <- ddf_satterthwaite_lambda(fit_phylolm_lambda, tree_rep)
  df_lambda <- res_satt_lambda$ddf
  dsatt_lambda <- c(dsatt_lambda, df_lambda)

  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
  pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = fit_phylolm$n - fit_phylolm$d))
  psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
  pbis <- c(pbis, 2 * pt(-abs(tval), df = ntips - 2))
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
### Satterthwaite BM REML
##############################################################################
##############################################################################
Nrep <- 100
pvanilla <- psatt <- psattapprox <- psatt_lambda <- pbis <- dsatt <- dsattapprox <- dsatt_lambda <- psatt_prod <- dsatt_prod <- NULL
effect <- 0
sigphylo <- 0.1
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  traits <- data.frame(species = species, id = species)

  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = sigphylo))

  traits <- data.frame(species = species,
                       id = species,
                       g1 = sim + rnorm(length(sim), 0, sd = sqrt(10)),
                       design = as.factor(design + 0), REML = TRUE)
  rownames(traits) <- species
  traits$g1[design] <- traits$g1[design] + effect

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM",
                                  measurement_error = TRUE,
                                  lower.bound = list(sigma2_error = getMinError(tree_rep)), REML = TRUE)
  # summary(fit_phylolm)
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
                                         measurement_error = FALSE,
                                         lower.bound = list(lambda = 1e-16),
                                         upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))), REML = TRUE)
  # fit_gls <- nlme::gls(g1 ~ design, traits, corPagel(0.5, tree_rep, form = ~species))

  res_satt <- ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)
  res_satt_approx <- ddf_satterthwaite_BM_error_approx(fit_phylolm, tree_rep)
  df <- res_satt$ddf[1,1]
  dsatt <- c(dsatt, df)
  dsattapprox <- c(dsattapprox, res_satt_approx$ddf)

  # res_satt_prod <- ddf_satterthwaite_BM(fit_phylolm, tree_rep, REML = TRUE)
  # df_prod <- res_satt_prod$df
  # dsatt_prod <- c(dsatt_prod, df_prod)
  #
  # res_satt_lambda <- ddf_satterthwaite_lambda(fit_phylolm_lambda, tree_rep, REML = TRUE)
  # df_lambda <- res_satt_lambda$ddf
  # dsatt_lambda <- c(dsatt_lambda, df_lambda)

  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]
  pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = fit_phylolm$n - fit_phylolm$d))
  psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
  psattapprox <- c(psattapprox, 2 * pt(-abs(tval), df = res_satt_approx$ddf))
  pbis <- c(pbis, 2 * pt(-abs(tval), df = ntips - 2))
  # psatt_lambda <- c(psatt_lambda, 2 * pt(-abs(tval), df = df_lambda))
  # psatt_prod <- c(psatt_prod, 2 * pt(-abs(tval), df = df_prod))
}
# hist(pvanilla)
# hist(psatt)
# hist(pbis)
# hist(psatt_lambda)
plot(sort(pvanilla))
points(sort(psatt), col = "red")
points(sort(psattapprox), col = "orange")
# points(sort(psatt_lambda), col = "orange")
points(sort(pbis), col = "green")
# points(sort(psatt_prod), col = "blue")
abline(a = 0, b = 1/Nrep)
plot(dsatt, ylim = c(0, fit_phylolm$n), col = "red")
points(dsattapprox, col = "orange")
# points(dsatt_lambda, col = "orange")
# points(dsatt_prod, col = "blue")
abline(a = fit_phylolm$n-2, b = 0)
abline(a = ntips-2, b = 0, col = "green")
mean(psatt <= 0.05)
mean(psattapprox <= 0.05)
# mean(psatt_lambda <= 0.05)
mean(pvanilla <= 0.05)
mean(pbis <= 0.05)

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
#   traits <- data.frame(species = species, id = species)
#
#   sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
#
#   traits <- data.frame(species = species,
#                        id = species,
#                        g1 = sim + rnorm(length(sim), 0, sd = sqrt(0.2)),
#                        design = as.factor(design + 0))
#   rownames(traits) <- species
#   traits$g1[design] <- traits$g1[design] + effect
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
#   approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/5))
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
#   Winv <- chol2inv(chol(W))
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
#   pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = n - d))
#   psatt <- c(psatt, 2 * pt(-abs(tval), df = df))
# }
# hist(pvanilla)
# hist(psatt)
# plot(sort(pvanilla))
# points(sort(psatt), col = "red")
# abline(a = 0, b = 1/Nrep)
# plot(dsatt)
# abline(a = n-2, b = 0)
# mean(psatt <= 0.05)
# mean(pvanilla <= 0.05)
#


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
#   traits <- data.frame(species = species, id = species)
#
#   sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
#
#   traits <- data.frame(species = species,
#                        id = species,
#                        g1 = sim + rnorm(length(sim), 0, sd = sqrt(0.2)),
#                        design = as.factor(design + 0))
#   rownames(traits) <- species
#   traits$g1[design] <- traits$g1[design] + effect
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
#     return(minusLogLik(x, y, X, tree_rep, "BM"))
#   }
#   approxHessian <- nlme::fdHess(pars = optpars, fun = fun, .relStep = .Machine$double.eps^(1/5))
#   A <- 1 / approxHessian$Hessian[1, 1] / (fit_phylolm$sigma2 / fit_phylolm$sigma2_error)^2
#
#   # approxHessian <- pracma::hessian(f = minusLogLik, x0 = fit_phylolm$sigma2_error / fit_phylolm$sigma2,
#   #                                  y = y, X = X, phy = tree_rep, model = "BM")
#   # A <- 1 / approxHessian[1, 1]
#
#   ## Satterthwaite
#   n <- length(tree_rep$tip.label)
#   d <- 2
#   K <- vcv(tree_rep)
#   Kd <- diag(diag(K))
#   V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
#   Vinv <- chol2inv(chol(V))
#   gamma <- fit_phylolm$sigma2_error / fit_phylolm$sigma2
#   W <- K + gamma * Kd
#   Winv <- chol2inv(chol(W))
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
#   derCgamma <- facmat %*% Kd %*% t(facmat)
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
# plot(dsatt)
# abline(a = n-2, b = 0)
# mean(psatt <= 0.05)
# mean(pvanilla <= 0.05)
#
