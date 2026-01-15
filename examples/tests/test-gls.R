# test_that("phylolm vs gls", {
#   set.seed(12891026)
#   ntips <- 10
#   tree <- ape::rphylo(ntips, 0.1, 0)
#   tree$edge.length <- tree$edge.length / vcv(tree)[1, 1]
#
#   r <- 3
#   traits <- data.frame(species = rep(tree$tip.label, r))
#   traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
#   rownames(traits) <- traits$id
#
#   tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)
#
#   traits <- traits[match(tree_rep$tip.label, traits$id), ]
#
#   design <- paste0("t", c(3, 8))
#   design <- design[!is.na(design)]
#   design <- traits$species %in% design
#   traits$design <- as.factor(design + 0)
#
#   plot(tree_rep)
#   tiplabels(pch = 21, col = traits$design, bg = traits$design)
#
#   sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.8))
#   traits$g1 <- sim + rnorm(length(sim), 0, sd = sqrt(0.2))
#   traits$g1[traits$design] <- traits$g1[traits$design] + effect
#
#   ###################
#   ## lambda
#   ###################
#
#   ## phylolm
#   fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
#                                          measurement_error = FALSE,
#                                          lower.bound = list(lambda = 1e-16),
#                                          upper.bound = list(lambda = getMaxLambda(getMinError(tree_rep))))
#   ## GLS
#   fit_gls <- nlme::gls(g1 ~ design, traits, correlation = corPagel(1, tree_rep, form = ~id), method = "ML")
#   ## phylolm bis
#   fit_phylolm_lambda <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "lambda",
#                                          measurement_error = FALSE,
#                                          lower.bound = list(lambda = fit_gls$modelStruct$corStruct[1]),
#                                          upper.bound = list(lambda = fit_gls$modelStruct$corStruct[1]),
#                                          starting.value = list(lambda = fit_gls$modelStruct$corStruct[1]))
#
#   expect_equal(fit_phylolm_lambda$coefficients, fit_gls$coefficients)
#   expect_equal(fit_phylolm_lambda$optpar, fit_gls$modelStruct$corStruct[1])
#   expect_equal(fit_phylolm_lambda$sigma2, fit_gls$sigma^2)
#
# })
