context("Rphylopars vs phylolm")

test_that("Rphylopars vs phylolm", {
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

  ## Check produced data
  rownames(traits) <- traits$id
  expect_warning(checkParamMatrix(traits, 'data', tree_rep, transpose = TRUE),
                 "`data` was not sorted in the correct order, when compared with the tips label. I am re-ordering it.")

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
