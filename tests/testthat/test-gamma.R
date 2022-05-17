##########################
# mass-extinction events #
##########################

test_that("Setting functions: constant mass-extinction event at time t.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setGammaConstant( 1, 0.5 )
  )

  # this should succeed
  expect_silent(
    tp$setGammaConstant( c(0,1), c(0.5, 0.5) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setGammaConstant( -1, 0.5 )
  )

  # this should fail (negative)
  expect_error(
    tp$setGammaConstant( 1, -0.5 )
  )

  # this should fail (too big)
  expect_error(
    tp$setGammaConstant( 1, 2.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaConstant( 1, c(0.5, 0.5) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaConstant( c(0, 1), c(0.5, 0.5, 0.5) )
  )

})


test_that("Setting functions: state-dependent mass-extinction event at time t", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setGammaStateDependent( c(1), matrix(c(0.5, 0.9), nrow = 1) )
  )

  # this should succeed
  expect_silent(
    tp$setGammaStateDependent( c(1, 2), matrix(c(0.5, 0.9, 0.1, 0.7), nrow = 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setGammaStateDependent( c(-1), matrix(c(0.5, 0.9), nrow = 1) )
  )

  # this should fail (negative)
  expect_error(
    tp$setGammaStateDependent( c(1), matrix(c(-0.5, 0.9), nrow = 1) )
  )

  # this should fail (too big)
  expect_error(
    tp$setGammaStateDependent( c(1), matrix(c(1.5, 0.9), nrow = 1) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaStateDependent( c(1), matrix(c(0.5, 0.9, 0.1, 0.7), nrow = 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaStateDependent( c(1), matrix(c(0.5, 0.9, 0.1), nrow = 1) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaStateDependent( c(1, 2), matrix(c(0.5, 0.9, 0.1, 0.7, 0.2, 0.8), nrow = 2) )
  )

})

test_that("Setting functions: constant mass-extinction event with state change at time t.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)
  P <- makeProbabilityMatrix(2)
  P[1,2] <- 0.5
  P[2,1] <- 0.2
  Ps <- c(P)
  PPs <- c(P, P)

  # this should succeed
  expect_silent(
    tp$setGammaAndZetaConstant(1, 0.5, Ps)
  )

  # this should succeed
  expect_silent(
    tp$setGammaAndZetaConstant( c(0,1), c(0.5, 0.5), PPs)
  )

  # this should fail (negative time)
  expect_error(
    tp$setGammaAndZetaConstant( -1, 0.5, Ps)
  )

  # this should fail (negative)
  expect_error(
    tp$setGammaAndZetaConstant( 1, -0.5, Ps)
  )

  # this should fail (too big)
  expect_error(
    tp$setGammaAndZetaConstant( 1, 2.0, Ps)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaAndZetaConstant( 1, c(0.5, 0.5), Ps)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaAndZetaConstant( c(0, 1), c(0.5, 0.5, 0.5), Ps)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setGammaAndZetaConstant( 1, 0.5, PPs)
  )

})
