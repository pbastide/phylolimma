context("addReplicates")

test_that("addReplicates", {
  skip_if_not_installed("phytools")
  ## Tree
  set.seed(1289)
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Simulate data
  reps <- sample(0:3, ntips, replace = TRUE)
  traits <- data.frame(g1 = rnorm(sum(reps), 1, 0.5),
                       species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps)))
  )
  traits$id <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  ## Replicates
  expect_warning(tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")
  expect_equal(length(tree_rep$tip.label), sum(reps))
  expect_equal(ape::is.ultrametric(tree_rep), TRUE)

  ## species number
  expect_equal(getSpeciesNumber(tree_rep), ntips - 2)

  ## Check produced data
  rownames(traits) <- traits$id
  expect_equal(traits, checkParamMatrix(traits, 'data', tree_rep, transpose = TRUE))

  ## Extra data
  traits_extra <- rbind(traits, data.frame(g1 = rnorm(3, 3, 3),
                                           species = rep("r1", 3),
                                           id = paste0("r1_", 1:3)))
  expect_error(addReplicatesOnTree(tree, traits_extra, species = "species", id = "id"),
               "are in the data but not in the tree.")

  ## Wrong names - species
  traits <- data.frame(g1 = rnorm(sum(reps), 1, 0.5),
                       labels = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps)))
  )
  traits$id <- mapply(paste0, traits$labels, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  expect_error(addReplicatesOnTree(tree, traits, species = "species", id = "id"),
               "data frame should contain a column named species with species names for each sample")
  expect_warning(addReplicatesOnTree(tree, traits, species = "labels", id = "id"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")

  ## Wrong names - ids
  traits <- data.frame(g1 = rnorm(sum(reps), 1, 0.5),
                       species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps)))
  )
  traits$unique_ids <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  expect_error(addReplicatesOnTree(tree, traits, species = "species", id = "id"),
               "data frame should contain a column named id with sample ids")
  expect_warning(addReplicatesOnTree(tree, traits, species = "species", id = "unique_ids"),
                 "Species 't4', 't9' are in the tree but not in the data. They will be droped from the final tree.")

  ## Replicated ids
  traits$unique_ids[2] <- traits$unique_ids[3]

  expect_error(addReplicatesOnTree(tree, traits, species = "species", id = "unique_ids"),
               "should be unique identifiers of the samples")
})

test_that("parseSpecies", {
  ## Tree
  set.seed(1289)
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Simulate data - missing species
  reps <- sample(0:3, ntips, replace = TRUE)
  traits <- data.frame(species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps))))
  traits$id <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  dat <- matrix(rnorm(sum(reps) * 10, 1, 0.5), nrow = 10)
  rownames(dat) <- paste0("g", 1:nrow(dat))
  colnames(dat) <- traits$id

  ids <- colnames(dat)
  expect_error(parse_species(tree, ids), "Sample ids could not be automatically matched")

  ## Simulate data - all good
  reps <- sample(1:3, ntips, replace = TRUE)
  traits <- data.frame(species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps))))
  traits$id <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  dat <- matrix(rnorm(sum(reps) * 10, 1, 0.5), nrow = 10)
  rownames(dat) <- paste0("g", 1:nrow(dat))
  colnames(dat) <- traits$id

  ids <- colnames(dat)
  expect_equal(parse_species(tree, ids),
               do.call(c, sapply(seq_along(tree$tip.label), function(i) rep(tree$tip.label[i], reps[i]))))

  ## Simulate data - wrong names
  reps <- sample(1:3, ntips, replace = TRUE)
  traits <- data.frame(species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps))))
  traits$id <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))

  dat <- matrix(rnorm(sum(reps) * 10, 1, 0.5), nrow = 10)
  rownames(dat) <- paste0("g", 1:nrow(dat))
  colnames(dat) <- traits$id
  colnames(dat) <- sub("t", "s", colnames(dat))

  ids <- colnames(dat)
  expect_error(parse_species(tree, ids), "Sample ids could not be automatically matched")

})

test_that("phylolmFit - rep", {
  skip_if_not_installed("phytools")

  ## Tree
  set.seed(1289)
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)

  ## Simulate data
  reps <- sample(1:5, ntips, replace = TRUE)
  traits <- data.frame(species = unname(do.call(c, mapply(function(x, y) rep(x, each = y), tree$tip.label, reps))))
  traits$id <- mapply(paste0, traits$species, do.call(c, sapply(reps, function(x) if (x > 0) paste0("_", 1:x))))
  dat <- matrix(rnorm(sum(reps) * 10, 1, 0.5), nrow = 10)
  rownames(dat) <- paste0("g", 1:nrow(dat))
  colnames(dat) <- traits$id
  nsamples <- ncol(dat)

  ## Design
  design <- matrix(1, nrow = nsamples, ncol = 2)
  colnames(design) <- c("(Intercept)", "condition")
  rownames(design) <- traits$id
  nullspecies <- paste0("t", sample(1:ntips, floor(ntips / 2)), "_")
  for (ns in nullspecies) {
    design[grep(ns, rownames(design)), 2] <- 0
  }

  ## Replicates
  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "id", eps = .Machine$double.eps^0.5)

  ## Fit with rep
  resPhyloLmFit <- phylolmFit(dat,
                              design = design,
                              phy = tree_rep,
                              model = "BM",
                              measurement_error = TRUE,
                              use_consensus = FALSE,
                              ddf_method = "Samples")

  ## Fit species tree - automatic
  resPhyloLmFitSpecies <- phylolmFit(dat,
                                     design = design,
                                     phy = tree,
                                     model = "BM",
                                     measurement_error = TRUE,
                                     use_consensus = FALSE,
                                     ddf_method = "Samples")

  expect_equal(resPhyloLmFit, resPhyloLmFitSpecies, tolerance = 1e-7)

  ## Fit species tree - col_species
  colnames(dat) <- sub("t", "s", colnames(dat))
  rownames(design) <- sub("t", "s", rownames(design))
  expect_error(phylolmFit(dat,
                          design = design,
                          phy = tree,
                          model = "BM",
                          measurement_error = TRUE,
                          use_consensus = FALSE,
                          ddf_method = "Samples"),
               "Sample ids could not be automatically matched against the species")

  resPhyloLmFitSpecies <- phylolmFit(dat,
                                     design = design,
                                     phy = tree,
                                     col_species = traits$species,
                                     model = "BM",
                                     measurement_error = TRUE,
                                     use_consensus = FALSE,
                                     ddf_method = "Samples")

  expect_equal(resPhyloLmFit$coefficients, resPhyloLmFitSpecies$coefficients, tolerance = 1e-7)
  expect_equal(resPhyloLmFit$sigma, resPhyloLmFitSpecies$sigma, tolerance = 1e-7)

})
