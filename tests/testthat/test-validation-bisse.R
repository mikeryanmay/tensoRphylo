#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: state-dependent birth-death model.", {

  # read tree data
  phy  <- readRDS(system.file("testdata", "extant_tree.Rda", package = "tensoRphylo"))
  data <- readRDS(system.file("testdata", "extant_data.Rda", package = "tensoRphylo"))

  # make parameters
  lambda <- c(0.1, 0.2)
  mu     <- c(0.03, 0.06)
  eta    <- 0.05

  # make a tp instance
  tp <- makeTensorPhylo(phy, data)

  # comparability: don't include the probability of the tree shape
  tp$setApplyTreeLikCorrection(FALSE)

  # set parameters
  tp$setLambdaStateDependent(lambda)
  tp$setMuStateDependent(mu)
  tp$setEtaConstantEqual(eta)

  # compute the likelihood
  ll <- tp$computeLogLikelihood()

  # compare tensorphylo against fixed value
  # (computed previously for the same dataset)
  expect_equal(ll, -36.1808415271)

  # make the diversitree model
  data_vec <- readRDS(system.file("testdata", "extant_data_vec.Rda", package = "tensoRphylo"))
  data_vec[1] <- NA
  mode(data_vec) <- "numeric"
  param_vec <- c(lambda, mu, eta, eta)
  bisse     <- diversitree::make.bisse(phy, data_vec)

  # compute the likelihood
  bisse_ll <- bisse(param_vec, root = 1, condition.surv = FALSE)

  # compare tensorphylo and diversitree
  expect_true(ll - bisse_ll < 1e-5)

})
