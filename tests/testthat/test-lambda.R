####################
# speciation rates #
####################

test_that("Setting functions: constant speciation rate.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setLambdaConstant( 1.0 )
  )

  # this should fail (negative)
  expect_error(
    tp$setLambdaConstant( -1.0 )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setLambdaConstant( c(1,2) )
  )

})

test_that("Setting functions: time-varying speciation rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setLambdaTimeDependent( c(1), c(1, 2) )
  )

  # this should fail (negative time)
  expect_error(
    tp$setLambdaTimeDependent( c(-1), c(1, 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setLambdaTimeDependent( c(1), c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setLambdaTimeDependent( c(1, 2), c(1, 2) )
  )

})

test_that("Setting functions: state-varying speciation rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setLambdaStateDependent( c(1,2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setLambdaStateDependent( c(-1, 2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setLambdaStateDependent( c(1, 2, 3) )
  )

})

test_that("Setting functions: state- and time-varying speciation rates.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setLambdaTimeStateDependent( c(1), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (negative rates)
  expect_error(
    tp$setLambdaTimeStateDependent( c(1), matrix(c(-1,2,-3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many times)
  expect_error(
    tp$setLambdaTimeStateDependent( c(1, 2), matrix(c(1,2,3,4), nrow = 2) )
  )

  # this should fail (dimensions: too many states)
  expect_error(
    tp$setLambdaTimeStateDependent( c(1), matrix(c(1,2,3,4,5,6), nrow = 2) )
  )

})
