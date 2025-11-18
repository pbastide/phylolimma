test_that("phylolm with transforms", {
  skip_if_not_installed("phylolm")
  skip_if_not_installed("phytools")

  ## Tree
  set.seed(1289)
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Replicates
  r <- 3
  ids <- as.vector(sapply(1:r, function(i) paste0(tree$tip.label, "_", i)))
  traits = data.frame(species = sub("\\_.", "", ids),
                      ids = ids)

  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "ids")

  ## traits
  sigma2_phylo <- 1
  sigma2_intra <- 0.1
  alpha <- 1
  resids <- phylolm::rTrait(n = 1,
                            phy = tree_rep,
                            model = "OU",
                            parameters = list(sigma2 = sigma2_phylo, alpha = alpha),
                            plot.tree = FALSE)
  resids <- resids + rnorm(ntips * r, mean = 0, sd = sqrt(sigma2_intra))
  traits <- data.frame(trait = resids)
  traits$species <- sub("\\_.", "", rownames(traits))
  traits <- traits[, c(2, 1)]
  traits$ids <- rownames(traits)

  ##############################################################################
  ### BM - lambda transform

  ## phylolm
  fit_phylolm <- phylolm::phylolm(trait~1, traits, tree_rep, measurement_error = TRUE)
  fit_phylolm_lambda <- phylolm::phylolm(trait~1, traits, tree_rep, model = "lambda",
                                         upper.bound = list(lambda = 1 - .Machine$double.eps))

  ## Test that lambda is the same as measurement error
  lambda_error <- fit_phylolm$sigma2 / (fit_phylolm$sigma2_error / max(ape::vcv(tree_rep)) + fit_phylolm$sigma2)
  expect_equal(lambda_error, fit_phylolm_lambda$optpar, tolerance = 1e-4)
  ## Likelihood
  expect_equivalent(fit_phylolm_lambda$logLik, fit_phylolm$logLik, tolerance = 1e-5)

  ## Scaled tree
  tree_bis <- rescale_tree(tree_rep)
  fit_phylolm_lambda_scaled <- phylolm::phylolm(trait~1, traits, tree_bis, model = "lambda",
                                                upper.bound = list(lambda = 1 - .Machine$double.eps))
  expect_equal(fit_phylolm_lambda_scaled$optpar, fit_phylolm_lambda$optpar, tolerance = 1e-5)
  expect_equivalent(fit_phylolm_lambda$logLik, fit_phylolm_lambda_scaled$logLik, tolerance = 1e-5)
  expect_equivalent(fit_phylolm_lambda$sigma2 * tree_height(tree_rep), fit_phylolm_lambda_scaled$sigma2, tolerance = 1e-5)


  ##############################################################################
  ### OU - lambda transform

  ## OU
  fit_ou <- phylolm::phylolm(trait ~ 1, traits, tree_rep, model = "OUfixedRoot", measurement_error = TRUE)

  ## Transform the tree
  tree_ou <- phylolm::transf.branch.lengths(tree_rep, "OUfixedRoot", parameters = list(alpha = fit_ou$optpar))
  tree_ou <- tree_ou$tree

  ## Lambda on the transformed tree
  fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_ou, model = "lambda", measurement_error = FALSE,
                                 upper.bound = list(lambda = 1 - .Machine$double.eps^0.5))

  ## lambda on OU
  tilde_t <- max(ape::vcv(tree_ou)) / (2 * fit_ou$optpar)
  lambda_ou_error <- fit_ou$sigma2 * tilde_t / (fit_ou$sigma2_error + fit_ou$sigma2 * tilde_t)

  ## Both are equal
  expect_equivalent(lambda_ou_error, fit_lambda$optpar, tolerance = 1e-5)

  ## Variances
  expect_equivalent(fit_ou$sigma2 / (2 * fit_ou$optpar) + fit_ou$sigma2_error / max(ape::vcv(tree_ou)), fit_lambda$sigma2, tolerance = 1e-4)

  ## Transform again
  tree_ou_lambda <- phylolm::transf.branch.lengths(tree_ou, "lambda", parameters = list(lambda = lambda_ou_error))
  tree_ou_lambda <- tree_ou_lambda$tree

  ## Fit BM on transformed tree
  fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_ou_lambda, model = "BM", measurement_error = FALSE)

  ## Same likelihood
  expect_equal(fit_bm$logLik, fit_ou$logLik)

  ## Variances
  expect_equivalent(fit_bm$sigma2, fit_lambda$sigma2, tolerance = 1e-4)

  ## Transform again
  tree_ou_error <- phylolm::transf.branch.lengths(tree_rep, "OUfixedRoot",
                                                  parameters = list(alpha = fit_ou$optpar,
                                                                    sigma2_error = 2 * fit_ou$optpar * fit_ou$sigma2_error / fit_ou$sigma2))
  tree_ou_error <- tree_ou_error$tree

  ## Fit BM on transformed tree
  fit_bm2 <- phylolm::phylolm(trait ~ 1, traits, tree_ou_error, model = "BM", measurement_error = FALSE)

  ## Same likelihood
  expect_equal(fit_bm2$logLik, fit_bm$logLik)

  ## Variance
  expect_equal(fit_bm2$sigma2 * max(diag(ape::vcv(tree_ou_error))),
               fit_bm$sigma2 * max(diag(ape::vcv(tree_ou_lambda))))

  ##############################################################################
  ### OU random root- lambda transform

  ## OU
  fit_ou <- phylolm::phylolm(trait ~ 1, traits, tree_rep, model = "OUrandomRoot", measurement_error = TRUE)

  ## Transform the tree
  tree_ou <- phylolm::transf.branch.lengths(tree_rep, "OUrandomRoot", parameters = list(alpha = fit_ou$optpar))
  tree_ou <- tree_ou$tree
  tree_ou$root.edge <- 0

  ## Lambda on the transformed tree
  fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_ou, model = "lambda", measurement_error = FALSE,
                                 upper.bound = list(lambda = 1 - .Machine$double.eps^0.5))

  ## lambda on OU
  tilde_t <- max(ape::vcv(tree_ou)) / (2 * fit_ou$optpar)
  lambda_ou_error <- fit_ou$sigma2 * tilde_t / (fit_ou$sigma2_error + fit_ou$sigma2 * tilde_t)

  ## Both are equal
  expect_equivalent(lambda_ou_error, fit_lambda$optpar, tolerance = 1e-5)

  ## Variances
  expect_equivalent(fit_ou$sigma2 / (2 * fit_ou$optpar) + fit_ou$sigma2_error / max(ape::vcv(tree_ou)), fit_lambda$sigma2, tolerance = 1e-4)

  ## Transform again
  tree_ou_lambda <- phylolm::transf.branch.lengths(tree_ou, "lambda", parameters = list(lambda = lambda_ou_error))
  tree_ou_lambda <- tree_ou_lambda$tree

  ## Fit BM on transformed tree
  fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_ou_lambda, model = "BM", measurement_error = FALSE)

  ## Same likelihood
  expect_equal(fit_bm$logLik, fit_ou$logLik)

  ## Variances
  expect_equivalent(fit_bm$sigma2, fit_lambda$sigma2, tolerance = 1e-4)

  ## Transform again
  tree_ou_error <- phylolm::transf.branch.lengths(tree_rep, "OUrandomRoot",
                                                  parameters = list(alpha = fit_ou$optpar,
                                                                    sigma2_error = 2 * fit_ou$optpar * fit_ou$sigma2_error / fit_ou$sigma2))
  tree_ou_error <- tree_ou_error$tree
  tree_ou_error$root.edge <- 0

  ## Fit BM on transformed tree
  fit_bm2 <- phylolm::phylolm(trait ~ 1, traits, tree_ou_error, model = "BM", measurement_error = FALSE)

  ## Same likelihood
  expect_equal(fit_bm2$logLik, fit_bm$logLik)

  ## Variance
  expect_equal(fit_bm2$sigma2 * max(diag(ape::vcv(tree_ou_error))),
               fit_bm$sigma2 * max(diag(ape::vcv(tree_ou_lambda))))


  ##############################################################################
  ### delta - lambda transform

  ## delta
  fit_delta <- phylolm::phylolm(trait ~ 1, traits, tree_rep,
                                model = "delta", upper.bound = list(delta = 100),
                                measurement_error = TRUE)

  ## Transform the tree
  tree_delta <- phylolm::transf.branch.lengths(tree_rep, "delta", parameters = list(delta = fit_delta$optpar))
  tree_delta <- tree_delta$tree

  ## Lambda on the transformed tree
  fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_delta, model = "lambda", measurement_error = FALSE)

  ## lambda on delta
  tilde_t <- max(ape::vcv(tree_delta))
  lambda_delta_error <- fit_delta$sigma2 * tilde_t / (fit_delta$sigma2_error + fit_delta$sigma2 * tilde_t)

  ## Both are equal
  expect_equivalent(lambda_delta_error, fit_lambda$optpar, tolerance = 1e-4)

  ## Transform again
  tree_delta_lambda <- phylolm::transf.branch.lengths(tree_delta, "lambda", parameters = list(lambda = lambda_delta_error))
  tree_delta_lambda <- tree_delta_lambda$tree

  ## Fit BM on transformed tree
  fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_delta_lambda, model = "BM", measurement_error = FALSE)

  ## Same likelihood
  expect_equal(fit_bm$logLik, fit_delta$logLik)

  ##############################################################################
  ### EB - lambda transform

  ## EB
  fit_EB <- phylolm::phylolm(trait ~ 1, traits, tree_rep, model = "EB", measurement_error = TRUE)

  ## Transform the tree
  tree_EB <- phylolm::transf.branch.lengths(tree_rep, "EB", parameters = list(EB = fit_EB$optpar))
  tree_EB <- tree_EB$tree

  ## Lambda on the transformed tree
  fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_EB, model = "lambda", measurement_error = FALSE,
                                 upper.bound = list(lambda = 1 - .Machine$double.eps))

  ## lambda on EB
  tilde_t <- max(ape::vcv(tree_EB))
  lambda_EB_error <- fit_EB$sigma2 * tilde_t / (fit_EB$sigma2_error + fit_EB$sigma2 * tilde_t)

  ## Both are equal
  expect_equivalent(lambda_EB_error, fit_lambda$optpar, tolerance = 1e-4)

  ## Transform again
  tree_EB_lambda <- phylolm::transf.branch.lengths(tree_EB, "lambda", parameters = list(lambda = lambda_EB_error))
  tree_EB_lambda <- tree_EB_lambda$tree

  ## Fit BM on transformed tree
  fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_EB_lambda, model = "BM", measurement_error = FALSE)

  ## Same likelihood
  expect_equal(fit_bm$logLik, fit_EB$logLik)

})

test_that("phylolm p-values", {
  skip_if_not_installed("Rphylopars")
  skip_if_not_installed("phylolm")

  ## Tree
  set.seed(1289)
  ntips <- 30
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## traits
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
  traits$x <- sample(c(0, 1), ntips * r, replace = TRUE)

  ## Replicates
  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id")

  ## phylolm
  fit_phylolm_error <- phylolm::phylolm(g1 ~ x, traits, tree_rep,
                                        measurement_error = TRUE, lower.bound = list(sigma2_error = .Machine$double.eps))
  fit_phylolm_lambda <- phylolm::phylolm(g1 ~ x, traits, tree_rep,
                                         model = "lambda", upper.bound = list(lambda = 1 - .Machine$double.eps),
                                         measurement_error = FALSE)

  expect_equal(summary(fit_phylolm_error)$coefficients[, "p.value"],
               summary(fit_phylolm_lambda)$coefficients[, "p.value"],
               tolerance = 1e-3)

})
