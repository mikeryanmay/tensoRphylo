##################################
# cladogenetic transition events #
##################################

test_that("Constructor functions: cladogenetic array", {

  # create the cladogenetic array
  Omega <- makeCladogeneticEvents(2)

  # this should succeed
  expect_silent(
    Omega[1,1,2] <- 0.1
  )

  # this should fail (negative)
  expect_error(
    Omega[1,1,2] <- -0.1
  )

  # this should fail (too big)
  expect_error(
    Omega[1,1,2] <- 1.2
  )

})

test_that("Setting functions: constant cladogenetic matrix.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)
  tp$setDebugMode(0)
  Omega <- makeCladogeneticEvents(num_states = 2)
  Omega[1,2,1] <- 0.25
  Omega[1,1,2] <- 0.25
  Omegas <- c(Omega, Omega)

  # this should succeed
  expect_silent(
    tp$setOmegaConstant(Omega)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setOmegaConstant( 1, Omegas  )
  )

})

test_that("Setting functions: time-varying cladogenetic matrix.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)
  tp$setDebugMode(0)
  Omega <- makeCladogeneticEvents(num_states = 2)
  Omega[1,2,1] <- 0.25
  Omega[1,1,2] <- 0.25
  Omegas <- c(Omega, Omega)

  # this should succeed
  expect_silent(
    tp$setOmegaTimeVarying(1, Omegas)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setOmegaTimeVarying( c(), Omegas  )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setOmegaTimeVarying( c(1, 2), Omegas  )
  )

})
