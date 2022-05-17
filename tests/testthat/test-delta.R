##############################
# destructive-sampling rates #
##############################

test_that("Setting functions: constant destructive-sampling rate.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

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
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaTimeDependent( c(1), c(1, 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setDeltaTimeDependent( c(-1), c(1, 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setDeltaTimeDependent( c(1), c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setDeltaTimeDependent( c(1, 2), c(1, 2) )
  )

})

test_that("Setting functions: state-varying destructive-sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaStateDependent( c(1,2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setDeltaStateDependent( c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setDeltaStateDependent( c(1, 2, 3) )
  )

})

test_that("Setting functions: state- and time-varying destructive-sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setDeltaTimeStateDependent( c(1), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setDeltaTimeStateDependent( c(1), matrix(c(-1,2,-3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many times)
  expect_error(
    tp$setDeltaTimeStateDependent( c(1, 2), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many states)
  expect_error(
    tp$setDeltaTimeStateDependent( c(1), matrix(c(1,2,3,4,5,6), nrow = 2) )
  )

})
