#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: time-dependent birth-death model.", {

  # read tree data
  phy <- tensoRphylo:::extant_tree

  # make parameters
  lambda   <- c(0.1, 0.2)
  lambda_t <- 10
  mu       <- 0.03

  # make a tp instance
  tp <- makeTensorPhylo(phy, nstates = 2)

  # set parameters
  tp$setLambdaTimeVarying(lambda_t, lambda)
  tp$setMuConstant(mu)

  # comparability: don't include the probability of the tree shape
  tp$setApplyTreeLikCorrection(FALSE)

  # comparability: condition on root survival
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  # compare tensorphylo against fixed value
  # (computed previously for the same dataset)
  expect_equal(
    tp$computeLogLikelihood(),
    -66.3502083005
  )

  # compare tensorphylo and TESS
  # (time is backward in TESS)
  bt <- branching.times(phy)
  expect_equal(
    tp$computeLogLikelihood(),
    TESS::tess.likelihood.rateshift( bt, rev(lambda), mu, max(bt) - lambda_t)
  )

})
