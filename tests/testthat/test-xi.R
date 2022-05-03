####################################
# mass-destructive-sampling events #
####################################

test_that("Setting functions: constant mass-destructive-sampling event at time t.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setXiConstant( 1, 0.5 )
  )

  # this should succeed
  expect_silent(
    tp$setXiConstant( c(0,1), c(0.5, 0.5) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setXiConstant( -1, 0.5 )
  )

  # this should fail (negative)
  expect_error(
    tp$setXiConstant( 1, -0.5 )
  )

  # this should fail (too big)
  expect_error(
    tp$setXiConstant( 1, 2.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setXiConstant( 1, c(0.5, 0.5) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setXiConstant( c(0, 1), c(0.5, 0.5, 0.5) )
  )

})


test_that("Setting functions: state-dependent mass-destructive-sampling event at time t", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setXiStateVarying( c(1), matrix(c(0.5, 0.9), nrow = 1) )
  )

  # this should succeed
  expect_silent(
    tp$setXiStateVarying( c(1, 2), matrix(c(0.5, 0.9, 0.1, 0.7), nrow = 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setXiStateVarying( c(-1), matrix(c(0.5, 0.9), nrow = 1) )
  )

  # this should fail (negative)
  expect_error(
    tp$setXiStateVarying( c(1), matrix(c(-0.5, 0.9), nrow = 1) )
  )

  # this should fail (too big)
  expect_error(
    tp$setXiStateVarying( c(1), matrix(c(1.5, 0.9), nrow = 1) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setXiStateVarying( c(1), matrix(c(0.5, 0.9, 0.1, 0.7), nrow = 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setXiStateVarying( c(1), matrix(c(0.5, 0.9, 0.1), nrow = 1) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setXiStateVarying( c(1, 2), matrix(c(0.5, 0.9, 0.1, 0.7, 0.2, 0.8), nrow = 2) )
  )

})
