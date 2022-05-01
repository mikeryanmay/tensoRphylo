##############################
# destructive-sampling rates #
##############################

test_that("Setting functions: constant destructive-sampling rate.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaConstant( 1.0 )
  )

  # this should fail (negative)
  expect_error(
    tp$setDeltaConstant( -1.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setDeltaConstant( c(1,2) )
  )

})

test_that("Setting functions: time-varying destructive-sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaTimeVarying( c(1), c(1, 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setDeltaTimeVarying( c(-1), c(1, 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setDeltaTimeVarying( c(1), c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setDeltaTimeVarying( c(1, 2), c(1, 2) )
  )

})

test_that("Setting functions: state-varying destructive-sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaStateVarying( c(1,2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setDeltaStateVarying( c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setDeltaStateVarying( c(1, 2, 3) )
  )

})

test_that("Setting functions: state- and time-varying destructive-sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaTimeStateVarying( c(1), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setDeltaTimeStateVarying( c(1), matrix(c(-1,2,-3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many times)
  expect_error(
    tp$setDeltaTimeStateVarying( c(1, 2), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many states)
  expect_error(
    tp$setDeltaTimeStateVarying( c(1), matrix(c(1,2,3,4,5,6), nrow = 2) )
  )

})










