##########################
# mass-speciation events #
##########################

test_that("Setting functions: constant mass-speciation event at time t.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setUpsilonConstant( 1, 0.5 )
  )

  # this should succeed
  expect_silent(
    tp$setUpsilonConstant( c(0,1), c(0.5, 0.5) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setUpsilonConstant( -1, 0.5 )
  )

  # this should fail (negative)
  expect_error(
    tp$setUpsilonConstant( 1, -0.5 )
  )

  # this should fail (too big)
  expect_error(
    tp$setUpsilonConstant( 1, 2.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setUpsilonConstant( 1, c(0.5, 0.5) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setUpsilonConstant( c(0, 1), c(0.5, 0.5, 0.5) )
  )

})


test_that("Setting functions: state-dependent mass-speciation event at time t", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setUpsilonStateVarying( c(1), matrix(c(0.5, 0.9), nrow = 1) )
  )

  # this should succeed
  expect_silent(
    tp$setUpsilonStateVarying( c(1, 2), matrix(c(0.5, 0.9, 0.1, 0.7), nrow = 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setUpsilonStateVarying( c(-1), matrix(c(0.5, 0.9), nrow = 1) )
  )

  # this should fail (negative)
  expect_error(
    tp$setUpsilonStateVarying( c(1), matrix(c(-0.5, 0.9), nrow = 1) )
  )

  # this should fail (too big)
  expect_error(
    tp$setUpsilonStateVarying( c(1), matrix(c(1.5, 0.9), nrow = 1) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setUpsilonStateVarying( c(1), matrix(c(0.5, 0.9, 0.1, 0.7), nrow = 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setUpsilonStateVarying( c(1), matrix(c(0.5, 0.9, 0.1), nrow = 1) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setUpsilonStateVarying( c(1, 2), matrix(c(0.5, 0.9, 0.1, 0.7, 0.2, 0.8), nrow = 2) )
  )

})
