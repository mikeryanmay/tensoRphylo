####################
# root frequencies #
####################

test_that("Setting functions: root frequencies.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # a functional example
  ## this should succeed
  expect_silent(
    tp$setRootPrior( 1:2 / 3  )
  )

  # dimensionality
  ## this should fail
  expect_error(
    tp$setRootPrior( rep(1, 3) / 3  )
  )

  # validity
  ## this should fail
  expect_error(
    tp$setRootPrior( c(1, 2) )
  )

})
