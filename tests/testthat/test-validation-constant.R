#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: constant-rate birth-death model.", {

  # read tree data
  phy <- tensoRphylo:::extant_tree

  # make parameters
  lambda <- 0.1
  mu     <- 0.03

  # make a tp instance
  tp <- makeTensorPhylo(phy, nstates = 2)

  # set parameters
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)

  # comparability: don't include the probability of the tree shape
  tp$setApplyTreeLikCorrection(FALSE)

  # comparability: condition on root survival
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  # compare tensorphylo against fixed value
  expect_equal(
    tp$computeLogLikelihood(),
    -66.8804540828
  )

  # compare tensorphylo and treepar
  expect_equal(
    tp$computeLogLikelihood(),
    -as.numeric(TreePar::LikConstant(lambda, mu, 1, TreeSim::getx(phy), survival = 1))
  )

  # compare tensorphylo and tess
  expect_equal(
    tp$computeLogLikelihood(),
    TESS::tess.likelihood(branching.times(phy), lambda = lambda, mu = mu)
  )

})
