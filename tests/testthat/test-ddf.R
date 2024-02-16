test_that("Satterthwaite - star tree", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")

  set.seed(1289)
  ## Tree
  ntips <- 10
  r <- 3
  tree <- ape::stree(ntips)
  tree$edge.length <- rep(1, nrow(tree$edge))

  traits <- data.frame(species = rep(tree$tip.label, r))
  traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
  rownames(traits) <- traits$id
  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)
  traits <- traits[match(tree_rep$tip.label, traits$id), ]

  ## Design
  design <- paste0("t", c(1))
  traits$design <- traits$species %in% design

  # plot(tree_rep)
  # ape::tiplabels(pch = 21, col = as.factor(traits$design + 0), bg = as.factor(traits$design + 0))

  ## Traits
  traits$g1 <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 10))
  traits$g1 <- traits$g1 + rnorm(length(traits$g1), mean = 0, sd = sqrt(0.2))
  traits$g1[traits$design] <- traits$g1[traits$design] + 0

  ## lmer reg
  fit_lmer <- lmerTest::lmer(g1 ~ 1 + design + (1|species), data = traits, REML = FALSE)

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE)

  ## Same estimates
  expect_equal(coef(summary(fit_lmer))[, 1], coef(summary(fit_phylolm))[, 1])
  expect_equal(lme4::VarCorr(fit_lmer)$species[1], fit_phylolm$sigma2, tolerance = 1e-7)
  expect_equal(fit_lmer@sigma^2, fit_phylolm$sigma2_error, tolerance = 1e-7)

  ## Satterthwaite
  expect_equal(anova(fit_lmer)$DenDF, ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)$ddf[1], tolerance = 1e-6)

  ## REML
  fit_lmer <- lmerTest::lmer(g1 ~ 1 + design + (1|species), data = traits, REML = TRUE)
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM", measurement_error = TRUE, REML = TRUE)
  expect_equal(coef(summary(fit_lmer))[, 1], coef(summary(fit_phylolm))[, 1])
  expect_equal(lme4::VarCorr(fit_lmer)$species[1], fit_phylolm$sigma2, tolerance = 1e-6)
  expect_equal(fit_lmer@sigma^2, fit_phylolm$sigma2_error, tolerance = 1e-7)
  expect_equal(anova(fit_lmer)$DenDF, ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)$ddf[1], tolerance = 1e-6)

})

test_that("Satterthwaite - chen tree", {
  ## Tree Chen small
  tree_rep <- ape::read.tree(text = "(((opossum.2:1.359619307e-14,opossum.1:0):1.359619307e-14,opossum.0:0):1,(((cow.2:1.359619307e-14,cow.1:0):1.359619307e-14,cow.0:0):0.6539768865,(((rat.2:1.359619307e-14,rat.1:0):1.359619307e-14,rat.0:0):0.6247450714,(human.1:1.359619307e-14,human.0:0):0.6247450714):0.02923181509):0.3460231135);")
  tree_rep$edge.length <- tree_rep$edge.length / ape::vcv(tree_rep)[1,1]
  species <- tree_rep$tip.label
  design <- paste0("rat.", 0:2)
  design <- species %in% design
  # plot(tree_rep)
  # ape::tiplabels(pch = 21, col = as.factor(design + 0), bg = as.factor(design + 0))
  ntips <- 4
  nsamples <- length(tree_rep$tip.label)

  ############################################################################
  ## Simulate data - High phylogenetic signal
  set.seed(1989)
  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 100))
  traits <- data.frame(species = species,
                       id = species,
                       g1 = sim + rnorm(length(sim), 0, sd = sqrt(0.1)),
                       design = as.factor(design + 0))
  rownames(traits) <- species
  traits$g1[design] <- traits$g1[design] + 0

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM",
                                  measurement_error = TRUE,
                                  lower.bound = list(sigma2_error = getMinError(tree_rep)))
  expect_equal(ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)$ddf[1], ntips, tolerance = 1e-2)

  ## phylolm - REML
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM",
                                  measurement_error = TRUE,
                                  lower.bound = list(sigma2_error = getMinError(tree_rep)),
                                  REML = TRUE)
  expect_equal(ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)$ddf[1], ntips - 2, tolerance = 1e-3)

  ############################################################################
  ## Simulate data - Low phylogenetic signal
  set.seed(1989)
  sim <- phylolm::rTrait(n = 1, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.01))
  traits <- data.frame(species = species,
                       id = species,
                       g1 = sim + rnorm(length(sim), 0, sd = sqrt(0.215)),
                       design = as.factor(design + 0))
  rownames(traits) <- species
  traits$g1[design] <- traits$g1[design] + 0

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM",
                                  measurement_error = TRUE,
                                  lower.bound = list(sigma2_error = getMinError(tree_rep)))
  expect_equal(ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)$ddf[1], nsamples, tolerance = 1e-2)

  ## phylolm - REML
  fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, tree_rep, model = "BM",
                                  measurement_error = TRUE,
                                  lower.bound = list(sigma2_error = getMinError(tree_rep)),
                                  REML = TRUE)
  expect_equal(ddf_satterthwaite_BM_error(fit_phylolm, tree_rep)$ddf[1], nsamples - 2, tolerance = 1e-2)

})

test_that("Satterthwaite - chen tree", {
  ## Tree Chen small
  tree_rep <- ape::read.tree(text = "(((opossum.2:1.359619307e-14,opossum.1:0):1.359619307e-14,opossum.0:0):1,(((cow.2:1.359619307e-14,cow.1:0):1.359619307e-14,cow.0:0):0.6539768865,(((rat.2:1.359619307e-14,rat.1:0):1.359619307e-14,rat.0:0):0.6247450714,(human.1:1.359619307e-14,human.0:0):0.6247450714):0.02923181509):0.3460231135);")
  tree_rep$edge.length <- tree_rep$edge.length / ape::vcv(tree_rep)[1,1]
  species <- tree_rep$tip.label
  design <- paste0("rat.", 0:2)
  design <- species %in% design
  # plot(tree_rep)
  # ape::tiplabels(pch = 21, col = as.factor(design + 0), bg = as.factor(design + 0))
  ntips <- 4
  nsamples <- length(tree_rep$tip.label)
  ngenes <- 20

  traits <- data.frame(species = species,
                       id = species,
                       design = as.factor(design + 0))
  rownames(traits) <- species

  design_matrix <- model.matrix(id ~ design, data = traits)

  ############################################################################
  ## Simulate data - High phylogenetic signal
  sim <- phylolm::rTrait(n = ngenes, tree_rep, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 100))
  y_data <- sim + rnorm(length(sim), 0, sd = sqrt(0.1))
  rownames(y_data) <- species
  y_data <- t(y_data)

  #################################################################################################
  ## With measurement error - Satterthwaite

  measurement_error <- TRUE
  ## Fit Phylo
  resPhyloLmFit <- phylolmFit(y_data, design = design_matrix, phy = tree_rep,
                              model = "BM",
                              measurement_error = TRUE,
                              use_consensus = FALSE,
                              ddf_method = "Satterthwaite")
  ## Fit Phylo Consensus
  resPhyloLmFitCons <- phylolmFit(y_data, design = design_matrix, phy = tree_rep,
                                  model = "BM",
                                  measurement_error = TRUE,
                                  use_consensus = TRUE,
                                  ddf_method = "Satterthwaite")

  ## Test ddf
  expect_equal(resPhyloLmFit$df.residual,resPhyloLmFitCons$df.residual)

  ## Test names and dimensions
  expect_equal(colnames(resPhyloLmFit$coefficients), colnames(resPhyloLmFitCons$coefficients))
  expect_equal(dim(resPhyloLmFit$coefficients), dim(resPhyloLmFitCons$coefficients))
  expect_equal(length(resPhyloLmFit$sigma), length(resPhyloLmFitCons$sigma))
  expect_equal(colnames(resPhyloLmFit$stdev.unscaled), colnames(resPhyloLmFitCons$stdev.unscaled))
  expect_equal(dim(resPhyloLmFit$stdev.unscaled), dim(resPhyloLmFitCons$stdev.unscaled))

  ## ebayes
  fitphy <- eBayes(resPhyloLmFit)
  fitphycons <- eBayes(resPhyloLmFitCons)

  ## Test names and dimensions
  expect_equal(colnames(fitphy$p.value), colnames(fitphycons$p.value))
  expect_equal(dim(fitphy$p.value), dim(fitphycons$p.value))

})
