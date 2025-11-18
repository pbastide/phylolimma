library(lmerTest)
library(phylolm)

source("examples/test_satterthwaite_utils.R")

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

## lmer reg
fit_lmer <- lmer(g1 ~ 1 + design + (1|species), data = traits, REML = FALSE)

## phylolm
fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep,
                                model = "BM", measurement_error = TRUE)

## Same estimates
all.equal(coef(summary(fit_lmer))[, 1], coef(summary(fit_phylolm))[, 1])
all.equal(VarCorr(fit_lmer)$species[1], fit_phylolm$sigma2)
all.equal(attr(VarCorr(fit_lmer), "sc")^2, fit_phylolm$sigma2_error)

## Satterthwaite
anova_lmer <- anova(fit_lmer)
res_phylo <- ddf_satterthwaite_sum(fit_phylolm, tree_rep, FALSE)
res_phylo_approx <- ddf_satterthwaite_sum_approx(fit_phylolm, tree_rep, FALSE)
F_phylo <- fit_phylolm$coefficients[2]^2 / (fit_phylolm$vcov[2, 2] * (ntot - 1) / ntot)

## Same degrees of freedom and stat values
all.equal(anova_lmer$DenDF, as.numeric(res_phylo$ddf))
all.equal(anova_lmer$`F value`, F_phylo)

## REML: should work better for small sample sizes
fit_lmer <- lmer(g1 ~ 1 + design + (1|species), data = traits, REML = TRUE)
fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE, REML = TRUE)
all.equal(coef(summary(fit_lmer))[, 1], coef(summary(fit_phylolm))[, 1])
all.equal(VarCorr(fit_lmer)$species[1], fit_phylolm$sigma2)
all.equal(attr(VarCorr(fit_lmer), "sc")^2, fit_phylolm$sigma2_error)

anova_lmer <- anova(fit_lmer)
res_phylo <- ddf_satterthwaite_sum(fit_phylolm, tree_rep, TRUE)
res_phylo_approx <- ddf_satterthwaite_sum_approx(fit_phylolm, tree_rep, TRUE)
F_phylo <- fit_phylolm$coefficients[2]^2 / fit_phylolm$vcov[2, 2]
all.equal(anova_lmer$DenDF, as.numeric(res_phylo$ddf))
all.equal(anova_lmer$`F value`, as.numeric(F_phylo))
