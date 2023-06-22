library(ape)
library(phytools)
library(mvMORPH)
library(phylolm)
library(MASS)

####################################################################################################
## simulate
####################################################################################################
set.seed(1289)

## tree
n <- 50
tree <- rphylo(n, 0.1, 0)
tree$edge.length <- tree$edge.length / max(diag(vcv(tree)))
plot(tree)
edgelabels()

## regimes
tree <- paintSubTree(tree, node = tree$edge[43, 2], state = "1", anc.state = "0", stem = TRUE)
tree <- paintSubTree(tree, node = tree$edge[91, 2], state = "1", anc.state = "0", stem = TRUE)
tree <- paintSubTree(tree, node = tree$edge[8, 2], state = "1", anc.state = "0", stem = TRUE)
plot(tree)

## Incidence matrices
Tfull <- PhylogeneticEM::incidence.matrix(tree)
Tfull <- Tfull %*% diag(sqrt(tree$edge.length))

V_tree <- vcv(tree)
all.equal(Tfull %*% t(Tfull), V_tree, check.attributes = FALSE)

## sub tree incidence matrix
edges_1 <- tree$mapped.edge[, 2] > 0
Tsub <- Tfull[, edges_1]

data_reg <- as.data.frame(cbind(rep(1, n), Tsub))
colnames(data_reg) <- c("intercept", paste0("B", which(edges_1)))
rownames(data_reg) <- tree$tip.label
qq <- qr(data_reg)
data_reg_qr <- data_reg[, qq$pivot[seq(qq$rank)]]

tree_sub <- tree
tree_sub$edge.length[!edges_1] <- 0
V_sub_tree <- vcv(tree_sub)
all.equal(Tsub %*% t(Tsub), V_sub_tree, check.attributes = FALSE)

anc <- 0

####################################################################################################
## Type I
####################################################################################################

sigma <- list("0" = 0.1,
              "1" = 0.1)
sigma_err <- 0.01

nrep <- 500
LTest <- vector(mode = "list", length = nrep)
Fstat <- rep(NA, nrep)
LeveneStatF <- rep(NA, nrep)
LeveneStatT <- rep(NA, nrep)
set.seed(1289)
for (rr in 1:nrep) {
  ## trait
  nsim <- 1
  ntraits <- 1
  traits <- mvSIM(tree,
                  nsim = nsim,
                  model = "BMM",
                  param = list(ntraits = ntraits,
                               sigma = sigma,
                               theta = anc))

  ## outlier
  traits[which.max(abs(traits))] <- 10 * traits[which.max(abs(traits))]

  ## measurement error
  traits <- traits + rnorm(n, mean = 0, sd = sqrt(sigma_err))

  ####################################################################################################
  ## fit
  ####################################################################################################

  ## Fit with error but same variance on the whole tree
  fit_error <- phylolm(traits ~ 1, phy = tree, model = "BM", measurement_error = TRUE)
  fit_pagel <- phylolm(traits ~ 1, phy = tree, model = "lambda", measurement_error = FALSE)
  lambda <- fit_pagel$optpar
  s_err <- fit_error$sigma2_error

  ## Fit with one regime
  fit_one <- mvBM(tree, traits, error = t(rep(s_err, n)), model = "BM1", method = "inverse", diagnostic = FALSE,echo = FALSE)
  ## test that this gives approx the same thing
  all.equal(fit_error$logLik, fit_one$LogLik)

  ## Fit with different regimes
  fit_two <- mvBM(tree, traits, error = t(rep(s_err, n)), model = "BMM", method = "inverse", diagnostic = FALSE,echo = FALSE)

  ####################################################################################################
  ## LR test
  ####################################################################################################

  LTest[[rr]] <- LRT(fit_two, fit_one, echo = FALSE)


  ####################################################################################################
  ## F test
  ####################################################################################################

  ## Metrix error
  Vlambda <- lambda * V_tree + (1 - lambda) * diag(rep(1, n))
  Clambda <- chol(Vlambda)
  # all.equal(t(Clambda) %*% Clambda, Vlambda, check.attributes = FALSE)

  ## transform data and incidence matrix
  trait_trans <- backsolve(Clambda, traits, transpose = TRUE)
  Tfull_trans <- backsolve(Clambda, Tfull, transpose = TRUE)
  Tsub_trans <- backsolve(Clambda, Tsub, transpose = TRUE)

  # all.equal(t(Tfull) %*% solve(Vlambda) %*% Tfull, t(Tfull_trans) %*% Tfull_trans, check.attributes = FALSE)

  ## F stat
  X <- backsolve(Clambda, rep(1, n), transpose = TRUE)
  PX <- X %*% solve(t(X) %*% X) %*% t(X)
  PXT <- diag(rep(1, n)) - PX
  SX <- t(trait_trans) %*% PXT %*% trait_trans

  SXbis <- (fit_one$sigma + s_err) * n
  SXter <- fit_pagel$sigma2 * n
  all.equal(SX, SXbis, check.attributes = FALSE, tol = 1e-5)
  all.equal(as.vector(SX), SXter, check.attributes = FALSE)

  W <- cbind(X, Tsub_trans)
  PW <- W %*% ginv(t(W) %*% W) %*% t(W)
  PWT <- diag(rep(1, n)) - PW
  SW <- t(trait_trans) %*% PWT %*% trait_trans

  # fit using big model and fixed lambda
  data_reg_trait <- data_reg_qr
  data_reg_trait$trait <- as.vector(traits)
  fixed_lambda <- fit_pagel$optpar
  fit_pagel_big <- phylolm(trait ~ . - 1, data = data_reg_trait, phy = tree, model = "lambda", measurement_error = FALSE, starting.value = fixed_lambda, lower.bound = fixed_lambda, upper.bound = fixed_lambda)

  SWbis <- fit_pagel_big$sigma2 * n
  all.equal(as.vector(SW), SWbis, check.attributes = FALSE, tol = 1e-2)

  p <- qr(X)$rank
  r <- qr(W)$rank

  Fstat[rr] <- ((SX - SW) / (r - p)) / (SW / (n - r))

  ####################################################################################################
  ## Levene test
  ####################################################################################################

  # Levene score
  Z <- abs(fit_pagel$residuals)
  # Proj on orthogonal of regression
  PXTraw <- diag(rep(1, n)) - fit_error$X %*% (t(fit_error$X) %*% fit_error$X)^{-1} %*% t(fit_error$X)
  # variance of the residuals UNDER H_0
  Vr <- PXTraw %*% Vlambda %*% PXTraw
  Cr <- cov2cor(Vr)
  # Expectation of z under H_0
  Elev0 <- sqrt(2 / pi) * sqrt(diag(Vr))
  # variance of levene score under H_0
  # Vlev <- matrix(NA, n, n)
  # for (i in 1:n) {
  #   for (j in 1:n) {
  #     # if (i != j) {
  #     #   Vlev[i,j] <- 2 / pi * Vr[i,j] * asin(Cr[i,j]) + 2 / pi * sqrt(Vr[i,i] * Vr[j,j] / (1 - Cr[i,j]^2)) - Elev0[i] * Elev0[j]
  #     # } else {
  #     #   Vlev[i,j] <- Vr[i,j] - Elev0[i]^2
  #     # }
  #     Vlev[i,j] <- 2 / pi * (sqrt(1 - Cr[i,j]^2) + Cr[i,j] * asin(Cr[i,j]) - 1 ) * sqrt(Vr[i,i] * Vr[j,j])
  #   }
  # }
  Vlev <- 2 / pi * (sqrt(1 - Cr^2) + Cr * asin(Cr) - 1 ) * sqrt(diag(Vr) %*% t(diag(Vr)))
  # additional variance of the residuals under H_1
  VlambdaSubtree <- lambda * V_sub_tree + (1 - lambda) * diag(rep(1, n))
  Vr1 <- PXTraw %*% VlambdaSubtree %*% PXTraw
  # Additional expectation under H_1
  Elev1 <- sqrt(2 / pi) * sqrt(diag(Vr1))

  ## F test on Z
  Clev <- chol(Vlev)
  Z_trans <- backsolve(Clev, Z, transpose = TRUE)
  Elev0_trans <- backsolve(Clev, Elev0, transpose = TRUE)
  Elev1_trans <- backsolve(Clev, Elev1, transpose = TRUE)

  fitlev0 <- lm(Z_trans ~ -1 + Elev0_trans)
  fitlev1 <- lm(Z_trans ~ -1 + Elev0_trans + Elev1_trans)
  LevAnova <- anova(fitlev0, fitlev1)
  LeveneStatF[rr] <- ifelse(is.na(LevAnova$F[2]), 0.0, LevAnova$F[2])

  LeveneStatT[rr] <- ifelse(is.na(fitlev1$coefficients[2]), 0.0, summary(fitlev1)$coefficients[2, "t value"])

}

## LRT
pvalLRT <- sapply(LTest, function(x) x$pval)
hist(pvalLRT)
sum(pvalLRT <= 0.05) / nrep

## F test
hist(Fstat, freq = FALSE)
xx <- seq(0, 10, 0.1)
lines(xx, df(xx, r - p, n - r))

pval <- pf(Fstat, r - p, n - r, lower.tail = FALSE)
hist(pval)
sum(pval <= 0.05) / nrep

## Levene test
hist(LeveneStatF, freq = FALSE)
xx <- seq(0, 10, 0.1)
lines(xx, df(xx, 2 - 1, n - 2))

pvallev <- pf(LeveneStatF, 2 - 1, n - 2, lower.tail = FALSE)
hist(pvallev)
sum(pvallev <= 0.05) / nrep

hist(LeveneStatT, freq = FALSE)
xx <- seq(-10, 10, 0.1)
lines(xx, dt(xx, n - 2))

pvallev <- pt(abs(LeveneStatT), n - 2, lower.tail = FALSE) * 2
hist(pvallev)
sum(pvallev <= 0.05) / nrep

####################################################################################################
## Type II
####################################################################################################

sigma <- list("0" = 0.1,
              "1" = 0.5)
sigma_err <- 0.01

nrep <- 500
LTest <- vector(mode = "list", length = nrep)
Fstat <- rep(NA, nrep)
LeveneStatF <- rep(NA, nrep)
LeveneStatT <- rep(NA, nrep)
set.seed(1289)
for (rr in 1:nrep) {
  ## trait
  nsim <- 1
  ntraits <- 1
  traits <- mvSIM(tree,
                  nsim = nsim,
                  model = "BMM",
                  param = list(ntraits = ntraits,
                               sigma = sigma,
                               theta = anc))

  ## outlier
  traits[which.max(abs(traits))] <- 10 * traits[which.max(abs(traits))]

  ## measurement error
  traits <- traits + rnorm(n, mean = 0, sd = sqrt(sigma_err))

  ####################################################################################################
  ## fit
  ####################################################################################################

  ## Fit with error but same variance on the whole tree
  fit_error <- phylolm(traits ~ 1, phy = tree, model = "BM", measurement_error = TRUE)
  fit_pagel <- phylolm(traits ~ 1, phy = tree, model = "lambda", measurement_error = FALSE)
  lambda <- fit_pagel$optpar
  s_err <- fit_error$sigma2_error

  ## Fit with one regime
  fit_one <- mvBM(tree, traits, error = t(rep(s_err, n)), model = "BM1", method = "inverse", diagnostic = FALSE, echo = FALSE)
  ## test that this gives approx the same thing
  all.equal(fit_error$logLik, fit_one$LogLik)

  ## Fit with different regimes
  fit_two <- mvBM(tree, traits, error = t(rep(s_err, n)), model = "BMM", method = "inverse", diagnostic = FALSE, echo = FALSE)

  ####################################################################################################
  ## LR test
  ####################################################################################################

  LTest[[rr]] <- LRT(fit_two, fit_one, echo = FALSE)


  ####################################################################################################
  ## F test
  ####################################################################################################

  ## Metrix error
  Vlambda <- lambda * V_tree + (1 - lambda) * diag(rep(1, n))
  Clambda <- chol(Vlambda)
  # all.equal(t(Clambda) %*% Clambda, Vlambda, check.attributes = FALSE)

  ## transform data and incidence matrix
  trait_trans <- backsolve(Clambda, traits, transpose = TRUE)
  Tfull_trans <- backsolve(Clambda, Tfull, transpose = TRUE)
  Tsub_trans <- backsolve(Clambda, Tsub, transpose = TRUE)

  # all.equal(t(Tfull) %*% solve(Vlambda) %*% Tfull, t(Tfull_trans) %*% Tfull_trans, check.attributes = FALSE)

  ## F stat
  X <- backsolve(Clambda, rep(1, n), transpose = TRUE)
  PX <- X %*% solve(t(X) %*% X) %*% t(X)
  PXT <- diag(rep(1, n)) - PX
  SX <- t(trait_trans) %*% PXT %*% trait_trans

  W <- cbind(X, Tsub_trans)
  PW <- W %*% ginv(t(W) %*% W) %*% t(W)
  PWT <- diag(rep(1, n)) - PW
  SW <- t(trait_trans) %*% PWT %*% trait_trans

  p <- qr(X)$rank
  r <- qr(W)$rank

  Fstat[rr] <- ((SX - SW) / (r - p)) / (SW / (n - r))

  ####################################################################################################
  ## Levene test
  ####################################################################################################

  # Levene score
  Z <- abs(fit_pagel$residuals)
  # Proj on orthogonal of regression
  PXTraw <- diag(rep(1, n)) - fit_error$X %*% (t(fit_error$X) %*% fit_error$X)^{-1} %*% t(fit_error$X)
  # variance of the residuals UNDER H_0
  Vr <- PXTraw %*% Vlambda %*% PXTraw
  Cr <- cov2cor(Vr)
  # Expectation of z under H_0
  Elev0 <- sqrt(2 / pi) * sqrt(diag(Vr))
  # variance of levene score under H_0
  # Vlev <- matrix(NA, n, n)
  # for (i in 1:n) {
  #   for (j in 1:n) {
  #     # if (i != j) {
  #     #   Vlev[i,j] <- 2 / pi * Vr[i,j] * asin(Cr[i,j]) + 2 / pi * sqrt(Vr[i,i] * Vr[j,j] / (1 - Cr[i,j]^2)) - Elev0[i] * Elev0[j]
  #     # } else {
  #     #   Vlev[i,j] <- Vr[i,j] - Elev0[i]^2
  #     # }
  #     Vlev[i,j] <- 2 / pi * (sqrt(1 - Cr[i,j]^2) + Cr[i,j] * asin(Cr[i,j]) - 1 ) * sqrt(Vr[i,i] * Vr[j,j])
  #   }
  # }
  Vlev <- 2 / pi * (sqrt(1 - Cr^2) + Cr * asin(Cr) - 1 ) * sqrt(diag(Vr) %*% t(diag(Vr)))
  # additional variance of the residuals under H_1
  VlambdaSubtree <- lambda * V_sub_tree + (1 - lambda) * diag(rep(1, n))
  Vr1 <- PXTraw %*% VlambdaSubtree %*% PXTraw
  # Additional expectation under H_1
  Elev1 <- sqrt(2 / pi) * sqrt(diag(Vr1))

  ## F test on Z
  Clev <- chol(Vlev)
  Z_trans <- backsolve(Clev, Z, transpose = TRUE)
  Elev0_trans <- backsolve(Clev, Elev0, transpose = TRUE)
  Elev1_trans <- backsolve(Clev, Elev1, transpose = TRUE)

  fitlev0 <- lm(Z_trans ~ -1 + Elev0_trans)
  fitlev1 <- lm(Z_trans ~ -1 + Elev0_trans + Elev1_trans)
  LevAnova <- anova(fitlev0, fitlev1)
  LeveneStatF[rr] <- ifelse(is.na(LevAnova$F[2]), 0.0, LevAnova$F[2])

  LeveneStatT[rr] <- ifelse(is.na(fitlev1$coefficients[2]), 0.0, summary(fitlev1)$coefficients[2, "t value"])

}

## LRT
pvalLRT <- sapply(LTest, function(x) x$pval)
hist(pvalLRT)
1 - sum(pvalLRT > 0.05) / nrep

## F test
hist(Fstat, freq = FALSE)
xx <- seq(0, 10, 0.1)
lines(xx, df(xx, r - p, n - r))

pval <- pf(Fstat, r - p, n - r, lower.tail = FALSE)
hist(pval)
1 - sum(pval > 0.05) / nrep

## Levene test
hist(LeveneStatF, freq = FALSE)
xx <- seq(0, 10, 0.1)
lines(xx, df(xx, 2 - 1, n - 2))

pvallev <- pf(LeveneStatF, 2 - 1, n - 2, lower.tail = FALSE)
hist(pvallev)
1 - sum(pvallev > 0.05) / nrep

hist(LeveneStatT, freq = FALSE)
xx <- seq(0, 10, 0.1)
lines(xx, dt(xx, n - 2))

pvallev <- pt(abs(LeveneStatT), n - 2, lower.tail = FALSE) * 2
hist(pvallev)
1 - sum(pvallev > 0.05) / nrep
