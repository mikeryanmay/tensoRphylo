#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: constant-rate serial-sampled birth-death model.", {

  # read tree data
  phy <- readRDS(test_path("testdata/serial_sampled_tree.Rda"))

  # make parameters
  lambda <- 0.1
  mu     <- 0.03
  delta  <- 0.05

  # make a tp instance
  tp <- makeTensorPhylo(phy, nstates = 2)
  tp$setApplyTreeLikCorrection(FALSE)

  # set parameters
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setDeltaConstant(delta)

  # comparability: condition on time
  tp$setConditionalProbabilityType(conditionalProbability$TIME)

  # compare tensorphylo against true value
  expect_equal(tp$computeLogLikelihood(), -40.6731823249)

})
