context("Rphylopars vs phylolm")

test_that("Rphylopars vs phylolm", {
  skip_if_not_installed("Rphylopars")
  skip_if_not_installed("phytools")

  ## Tree
  set.seed(1289)
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Simulate data
  r <- 3
  traits <- Rphylopars::simtraits(tree = tree, ntraits = 1, nreps = r, nmissing = 0,
                                  v = 1, anc = 0,
                                  intraspecific = 0.1,
                                  model="delta", parameters = list(alpha = 1.5),
                                  nsim = 1)
  traits <- traits$trait_data
  colnames(traits)[2] <- "g1"
  traits$id <- mapply(paste0, traits$species, paste0("_", rep(1:r, each = ntips)))
  rownames(traits) <- traits$id

  ## Replicates
  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)

  ##############################################################################
  ### BM

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ 1, traits, tree_rep, measurement_error = TRUE)

  ## Rphylopars
  fit_phylopars <- Rphylopars::phylopars(traits[, 1:2], tree, REML = FALSE)

  ## Test phylo variance
  expect_equivalent(fit_phylopars$pars$phylocov, fit_phylolm$sigma2, tolerance = 1e-6)
  ## Test pheno variance
  expect_equivalent(fit_phylopars$pars$phenocov, fit_phylolm$sigma2_error, tolerance = 1e-7)
  ## Likelihood
  expect_equivalent(fit_phylopars$logLik, fit_phylolm$logLik)

  ## Timing
  # trait_phylopars <- traits[, 1:2]
  # microbenchmark::microbenchmark(phylolm::phylolm(g1 ~ 1, traits, tree_rep, measurement_error = TRUE),
  #                                Rphylopars::phylopars(trait_phylopars, tree, REML = FALSE))

  ##############################################################################
  ### OU

  ## phylolm
  fit_phylolm <- phylolm::phylolm(g1 ~ 1, traits, tree_rep,
                                  model = "OUfixedRoot", measurement_error = TRUE)

  ## Rphylopars
  fit_phylopars <- Rphylopars::phylopars(traits[, 1:2], tree, model = "OU", REML = FALSE)

  ## Test phylo variance
  expect_equivalent(fit_phylopars$pars$phylocov * 2 * fit_phylopars$model$alpha, fit_phylolm$sigma2, tolerance = 1e-3)
  ## Test pheno variance
  expect_equivalent(fit_phylopars$pars$phenocov, fit_phylolm$sigma2_error, tolerance = 1e-5)
  ## Test alpha
  expect_equivalent(fit_phylopars$model$alpha, fit_phylolm$optpar, tolerance = 1e-4)
  ## Likelihood
  expect_equivalent(fit_phylopars$logLik, fit_phylolm$logLik, tolerance = 1e-6)

  ## Timing
  # trait_phylopars <- traits[, 1:2]
  # microbenchmark::microbenchmark(phylolm::phylolm(g1 ~ 1, traits, tree_rep,
  #                                                 model = "OUfixedRoot", measurement_error = TRUE,
  #                                                 lower.bound = list(sigma2_error = 1e-7)),
  #                                Rphylopars::phylopars(trait_phylopars, tree, model = "OU", REML = FALSE))

  ##############################################################################
  ### With predictors - BM

  # ## Condition
  # traits$cond <- rep(sample(c(0, 1), ntips, replace = TRUE), r)
  # traits[traits$cond == 1, "g1"] <- traits[traits$cond == 1, "g1"] + 5 ## strong effect
  #
  # ## phylolm
  # fit_phylolm <- phylolm::phylolm(g1 ~ 1 + cond, traits, tree_rep, measurement_error = TRUE)
  #
  # ## Rphylopars
  # fit_phylopars <- Rphylopars::phylopars.lm(g1 ~ 1 + cond, traits[, c(1:2, 4)], tree, pheno_error = TRUE, REML = FALSE)
  #
  # ## Test phylo variance
  # expect_equivalent(fit_phylopars$pars$phylocov, fit_phylolm$sigma2, tolerance = 1e-6)
  # ## Test pheno variance
  # expect_equivalent(fit_phylopars$pars$phenocov, fit_phylolm$sigma2_error, tolerance = 1e-7)
  # ## Likelihood
  # expect_equivalent(fit_phylopars$logLik, fit_phylolm$logLik)

})



