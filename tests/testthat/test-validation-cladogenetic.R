#####################################
# validation (checking likelihoods) #
#####################################

test_that("Validation: cladogenetic model.", {

  # read tree data
  phy  <- readRDS(test_path("testdata/extant_tree.Rda"))
  data <- readRDS(test_path("testdata/extant_data.Rda"))

  # parameters
  lambda <- 1
  mu     <- 0.5
  q      <- 0.1
  pc     <- c(0.4, 0.3)
  pa     <- c(0.7, 0.2)

  # make a tp instance
  tp <- makeTensorPhylo(phy, data)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$TIME)

  # make the tp parameters
  Q <- makeRateMatrix(2, q)
  O <- makeCladogeneticEvents(2)
  O[1,1,2] <- pc[1] * pa[1] * 0.5
  O[1,2,1] <- pc[1] * pa[1] * 0.5
  O[1,2,2] <- pc[1] * (1 - pa[1])
  O[2,2,1] <- pc[2] * pa[2] * 0.5
  O[2,1,2] <- pc[2] * pa[2] * 0.5
  O[2,1,1] <- pc[2] * (1 - pa[2])

  # set the parameters
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setEtaConstantUnequal(Q)
  tp$setOmegaConstant(O)

  # switch to classe params
  pars <- diversitree::starting.point.classe(phy, 2)
  pars["lambda111"] <- lambda * (1 - pc[1])
  pars["lambda112"] <- lambda * pc[1] * pa[1]
  pars["lambda122"] <- lambda * pc[1] * (1 - pa[1])
  pars["lambda211"] <- lambda * pc[2] * (1 - pa[2])
  pars["lambda212"] <- lambda * pc[2] * pa[2]
  pars["lambda222"] <- lambda * (1 - pc[2])
  pars["mu1"] <- mu
  pars["mu2"] <- mu
  pars["q12"] <- q
  pars["q21"] <- q

  # make the classe model
  data_vec  <- (data %*% c(0,1))[,1] + 1
  model <- diversitree::make.classe(phy, data_vec, k=2, sampling.f=c(1,1))

  # compare tensorphylo against fixed value
  # (computed previously for the same dataset)
  expect_equal(
    tp$computeLogLikelihood(),
    -42.7394081639
  )

  # compare tensorphylo and diversitree
  # these won't be exactly the same because of numerical details
  expect_true(
    tp$computeLogLikelihood() - model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT) < 1e-4
  )

})
