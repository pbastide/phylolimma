# test_that("initialisation", {
#   skip_if_not_installed("phylolm")
#
#   set.seed(12891026)
#   ## Tree
#   ntips <- 100
#   tree <- ape::rphylo(ntips, 0.1, 0)
#   mat_tree <- ape::vcv(tree)
#   ## params
#   alpha_true <- 0.5
#   gamma_true <- 1 / 2 / 0.5
#   ## data
#   y_data <- phylolm::rTrait(1, tree, model = "OU", parameters = list(alpha = 0.5))
#   ## Design
#   design <- matrix(1, nrow = ntips, ncol = 2)
#   design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
#   colnames(design) <- c("(Intercept)", "condition")
#   rownames(design) <- tree$tip.label
#   y_data <- design %*% c(3, -2) + y_data
#
#   ## gamma parameter
#   gamma_init <- initGammaParameter(design, y_data, tree)
#   expect_equal(gamma_init, 0.9243532, tol = 1e-3)
#
#   ## alpha parameter
#
# })
