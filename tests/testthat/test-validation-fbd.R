#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: constant-rate fossilized birth-death model.", {

  # read tree data
  phy <- readRDS(system.file("testdata", "sampled_ancestor_tree.Rda", package = "tensoRphylo"))

  # make parameters
  lambda <- 0.1
  mu     <- 0.03
  phi    <- 0.05

  # make a tp instance
  tp <- makeTensorPhylo(phy, nstates = 2)

  # set parameters
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setPhiConstant(phi)

  # comparability: don't include the probability of the tree shape
  tp$setApplyTreeLikCorrection(FALSE)

  # compare tensorphylo against true value
  expect_equal(tp$computeLogLikelihood(), -58.1330970357)

})
