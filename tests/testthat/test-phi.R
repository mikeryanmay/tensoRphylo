####################
# sampling rates #
####################

test_that("Setting functions: constant sampling rate.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setPhiConstant( 1.0 )
  )

  # this should fail (negative)
  expect_error(
    tp$setPhiConstant( -1.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setPhiConstant( c(1,2) )
  )

})

test_that("Setting functions: time-varying sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setPhiTimeVarying( c(1), c(1, 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setPhiTimeVarying( c(-1), c(1, 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setPhiTimeVarying( c(1), c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setPhiTimeVarying( c(1, 2), c(1, 2) )
  )

})

test_that("Setting functions: state-varying sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setPhiStateVarying( c(1,2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setPhiStateVarying( c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setPhiStateVarying( c(1, 2, 3) )
  )

})

test_that("Setting functions: state- and time-varying sampling rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setPhiTimeStateVarying( c(1), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setPhiTimeStateVarying( c(1), matrix(c(-1,2,-3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many times)
  expect_error(
    tp$setPhiTimeStateVarying( c(1, 2), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many states)
  expect_error(
    tp$setPhiTimeStateVarying( c(1), matrix(c(1,2,3,4,5,6), nrow = 2) )
  )

})
