####################
# extinction rates #
####################

test_that("Setting functions: constant extinction rate.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setMuConstant( 1.0 )
  )

  # this should fail (negative)
  expect_error(
    tp$setMuConstant( -1.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setMuConstant( c(1,2) )
  )

})

test_that("Setting functions: time-varying extinction rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setMuTimeVarying( c(1), c(1, 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setMuTimeVarying( c(-1), c(1, 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setMuTimeVarying( c(1), c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setMuTimeVarying( c(1, 2), c(1, 2) )
  )

})

test_that("Setting functions: state-varying extinction rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setMuStateVarying( c(1,2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setMuStateVarying( c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setMuStateVarying( c(1, 2, 3) )
  )

})

test_that("Setting functions: state- and time-varying extinction rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhylo, 2)

  # this should succeed
  expect_silent(
    tp$setMuTimeStateVarying( c(1), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setMuTimeStateVarying( c(1), matrix(c(-1,2,-3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many times)
  expect_error(
    tp$setMuTimeStateVarying( c(1, 2), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many states)
  expect_error(
    tp$setMuTimeStateVarying( c(1), matrix(c(1,2,3,4,5,6), nrow = 2) )
  )

})










