##########################
# mass-extinction events #
##########################

test_that("Setting functions: constant mass-state-change event at time t.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)
  P <- makeProbabilityMatrix(2)
  P[1,2] <- 0.5
  P[2,1] <- 0.2
  Ps <- c(P)
  PPs <- c(P, P)

  # this should succeed
  expect_silent(
    tp$setZeta(1, c(P))
  )

  # this should succeed
  expect_silent(
    tp$setZeta( c(0,1), c(P,P))
  )

  # this should fail (negative time)
  expect_error(
    tp$setZeta( -1, c(P))
  )

  # this should fail (dimensions)
  expect_error(
    tp$setZeta( c(1,2), Ps)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setZeta( 1, PPs)
  )

})
