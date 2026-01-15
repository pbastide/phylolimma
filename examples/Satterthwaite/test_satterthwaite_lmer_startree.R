library(lmerTest)
library(phylolm)

set.seed(12891289)
## Tree
ntips <- 4
r <- 3
tree <- stree(ntips)
tree$edge.length <- rep(1, nrow(tree$edge))

traits <- data.frame(species = rep(tree$tip.label, r))
traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
rownames(traits) <- traits$id
tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)
traits <- traits[match(tree_rep$tip.label, traits$id), ]

design <- paste0("t", c(1))
traits$design <- traits$species %in% design

plot(tree_rep)
tiplabels(pch = 21, col = as.factor(traits$design + 0), bg = as.factor(traits$design + 0))

## Traits
traits$g1 <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
traits$g1 <- traits$g1 + rnorm(length(traits$g1), mean = 0, sd = sqrt(0.2))
traits$g1[traits$design] <- traits$g1[traits$design] + 0

## lmer reg
fit_lmer <- lmer(g1 ~ 1 + design + (1|species), data = traits, REML = FALSE)

## phylolm
fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE)
summary(fit_phylolm)

## Same estimates
all.equal(coef(summary(fit_lmer))[, 1], coef(summary(fit_phylolm))[, 1])
all.equal(VarCorr(fit_lmer)$species[1], fit_phylolm$sigma2)
all.equal(attr(VarCorr(fit_lmer), "sc")^2, fit_phylolm$sigma2_error)

## Satterthwaite
anova(fit_lmer)

ddf_satterthwaite_sum(fit_phylolm, tree_rep, FALSE)
fit_phylolm$coefficients[1] / sqrt(res_lm$vcov[1, 1])

ddf_satterthwaite_BM(fit_phylolm, tree_rep, FALSE)

## REML
## lmer reg
fit_lmer <- lmer(g1 ~ 1 + design + (1|species), data = traits, REML = TRUE)
fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE, REML = TRUE)
all.equal(coef(summary(fit_lmer))[, 1], coef(summary(fit_phylolm))[, 1])
all.equal(VarCorr(fit_lmer)$species[1], fit_phylolm$sigma2)
all.equal(attr(VarCorr(fit_lmer), "sc")^2, fit_phylolm$sigma2_error)
anova(fit_lmer)
ddf_satterthwaite_sum(fit_phylolm, tree_rep, REML = TRUE)
fit_phylolm$coefficients[2]^2 / fit_phylolm$vcov[2, 2]

ddf_satterthwaite_BM(fit_phylolm, tree_rep, TRUE)

##############################################################################
##############################################################################
### Null model
##############################################################################
##############################################################################
Nrep <- 100
pvanilla <- psattsum <- psattprod <- plmer <- pbis <- dsattsum <- dsattprod <- dlmer <- psatt_lambda <- dsatt_lambda <- NULL
effect <- 0
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  lam <- 0.5
  traits$g1 <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = lam))
  traits$g1 <- traits$g1 + rnorm(length(traits$g1), mean = 0, sd = sqrt(1-lam))
  traits$g1[traits$design] <- traits$g1[traits$design] + 0

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE, lower.bound = list(lambda = 1e-16))
  # summary(fit_phylolm)
  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]

  fit_lmer <- lmer(g1 ~ 1 + design + (1|species), data = traits, REML = FALSE)

  res_satt <- ddf_satterthwaite_sum(fit_phylolm, tree_rep)
  dfsattsum <- res_satt$ddf
  dsattsum <- c(dsattsum, dfsattsum)

  # res_satt <- ddf_satterthwaite_BM(fit_phylolm, tree_rep)
  # dfsattprod <- res_satt$df
  # dsattprod <- c(dsattprod, dfsattprod)

  ## phylolm lambda
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
                                         measurement_error = FALSE,
                                         lower.bound = list(lambda = 1e-16),
                                         upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))))

  res_satt_lambda <- ddf_satterthwaite_lambda(fit_phylolm_lambda, tree_rep)
  df_lambda <- res_satt_lambda$ddf
  dsatt_lambda <- c(dsatt_lambda, df_lambda)

  dflmer <- anova(fit_lmer)$DenDF
  dlmer <- c(dlmer, dflmer)

  pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = fit_phylolm$n - fit_phylolm$d))
  psattsum <- c(psattsum, 2 * pt(-abs(tval), df = dfsattsum))
  # psattprod <- c(psattprod, 2 * pt(-abs(tval), df = dfsattprod))
  plmer <- c(plmer, 2 * pt(-abs(tval), df = dflmer))
  pbis <- c(pbis, 2 * pt(-abs(tval), df = ntips - 2))
  psatt_lambda <- c(psatt_lambda, 2 * pt(-abs(tval), df = df_lambda))
}

all.equal(plmer, psattsum)

plot(sort(pvanilla))
points(sort(plmer), col = "blue")
points(sort(psattsum), col = "red")
points(sort(psatt_lambda), col = "orange")
points(sort(psattprod), col = "yellow")
points(sort(pbis), col = "green")
abline(a = 0, b = 1/Nrep)
plot(dlmer, col = "blue", ylim = c(0, 15))
points(dsattsum, col = "red")
points(dsatt_lambda, col = "orange")
points(dsattprod, col = "yellow")
abline(a = ntips * r - 2, b = 0)
abline(a = ntips - 2, b = 0, col = "green")
mean(pvanilla <= 0.05)
mean(psattsum <= 0.05)
mean(psattprod <= 0.05)
mean(psatt_lambda <= 0.05)
mean(plmer <= 0.05)
mean(pbis <= 0.05)

##############################################################################
##############################################################################
### Null model - REML
##############################################################################
##############################################################################
Nrep <- 100
pvanilla <- psattsum <- psattprod <- plmer <- pbis <- dsattsum <- dsattprod <- dlmer <- psatt_lambda <- dsatt_lambda <- NULL
effect <- 0
set.seed(1289)
for (rep in 1:Nrep) {

  ## Simulate data
  lam <- 0.5
  traits$g1 <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = lam))
  traits$g1 <- traits$g1 + rnorm(length(traits$g1), mean = 0, sd = sqrt(1-lam))
  traits$g1[traits$design] <- traits$g1[traits$design] + 0

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE, lower.bound = list(lambda = 1e-16), REML = TRUE)
  # summary(fit_phylolm)
  tval <- summary(fit_phylolm)$coefficients[2, "t.value"]

  fit_lmer <- lmer(g1 ~ 1 + design + (1|species), data = traits, REML = TRUE)

  res_satt <- ddf_satterthwaite_sum(fit_phylolm, tree_rep, REML = TRUE)
  dfsattsum <- res_satt$ddf
  dsattsum <- c(dsattsum, dfsattsum)

  # res_satt <- ddf_satterthwaite_BM(fit_phylolm, tree_rep)
  # dfsattprod <- res_satt$df
  # dsattprod <- c(dsattprod, dfsattprod)

  ## phylolm lambda
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
                                         measurement_error = FALSE,
                                         lower.bound = list(lambda = 1e-16),
                                         upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))),
                                         REML = TRUE)

  res_satt_lambda <- ddf_satterthwaite_lambda(fit_phylolm_lambda, tree_rep, REML = TRUE)
  df_lambda <- res_satt_lambda$ddf
  dsatt_lambda <- c(dsatt_lambda, df_lambda)

  dflmer <- anova(fit_lmer)$DenDF
  dlmer <- c(dlmer, dflmer)

  pvanilla <- c(pvanilla, 2 * pt(-abs(tval), df = fit_phylolm$n - fit_phylolm$d))
  psattsum <- c(psattsum, 2 * pt(-abs(tval), df = dfsattsum))
  # psattprod <- c(psattprod, 2 * pt(-abs(tval), df = dfsattprod))
  plmer <- c(plmer, 2 * pt(-abs(tval), df = dflmer))
  pbis <- c(pbis, 2 * pt(-abs(tval), df = ntips - 2))
  psatt_lambda <- c(psatt_lambda, 2 * pt(-abs(tval), df = df_lambda))
}

all.equal(dlmer, dsattsum)

plot(sort(pvanilla))
points(sort(plmer), col = "blue")
points(sort(psattsum), col = "red")
points(sort(psatt_lambda), col = "orange")
# points(sort(psattprod), col = "yellow")
points(sort(pbis), col = "green")
abline(a = 0, b = 1/Nrep)
plot(dlmer, col = "blue", ylim = c(0, 15))
points(dsattsum, col = "red")
points(dsatt_lambda, col = "orange")
# points(dsattprod, col = "yellow")
abline(a = ntips * r - 2, b = 0)
abline(a = ntips - 2, b = 0, col = "green")
mean(pvanilla <= 0.05)
mean(psattsum <= 0.05)
mean(psatt_lambda <= 0.05)
mean(plmer <= 0.05)
mean(pbis <= 0.05)
