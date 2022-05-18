#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: constant-rate birth-death model with character evolution.", {

  # read tree data
  phy  <- readRDS(system.file("testdata", "extant_tree.Rda", package = "tensoRphylo"))
  data <- readRDS(system.file("testdata", "extant_data.Rda", package = "tensoRphylo"))

  # don't want to deal with missing data for phytools
  data[1,1] <- 0

  # make parameters
  lambda <- 0.1
  mu     <- 0.03
  Q      <- makeRateMatrix(2)
  Q[1,2] <- 0.01
  Q[2,1] <- 0.02

  # make a tp instance
  tp <- makeTensorPhylo(phy, data)

  # set parameters
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setEtaConstantUnequal(Q)

  # comparability: don't include the probability of the tree shape
  tp$setApplyTreeLikCorrection(FALSE)

  # comparability: condition on root survival
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  # compare tensorphylo against fixed value
  # (computed previously for the same dataset)
  expect_equal(
    tp$computeLogLikelihood(),
    -38.6468916843
  )

  # compare tensorphylo and treepar plus phytools
  treepar_ll  <- -as.numeric(TreePar::LikConstant(lambda, mu, 1, TreeSim::getx(phy), survival = 1))
  phytools_ll <- as.numeric(phytools:::getPars(phy, data, "ARD", Q$getMatrix(), phy, 1e-16, 2, pi = c(0.5,0.5))$loglik)
  expect_equal(tp$computeLogLikelihood(), treepar_ll + phytools_ll)

})
