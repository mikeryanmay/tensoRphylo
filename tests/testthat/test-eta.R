################################
# anagenetic transition events #
################################

test_that("Constructor functions: anagenetic matrix", {

  # create the anagenetic matrix
  Q <- makeRateMatrix(2)

  # this should succeed
  expect_silent(
    Q[1,2] <- 0.1
  )

  # this should fail (negative)
  expect_error(
    Q[1,2] <- -0.1
  )

  # this should fail (dimensions)
  expect_error(
    Q[1,3] <- 0.1
  )

})

test_that("Setting functions: constant equal anagenetic matrix.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setEtaConstantEqual(0.1)
  )

  # this should fail (negative)
  expect_error(
    tp$setEtaConstantEqual(-0.1)
  )

  # this should fail (dimensions)
  expect_error(
    tp$setEtaConstantEqual(c(1,2))
  )

})

test_that("Setting functions: constant unequal anagenetic matrix.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)
  Q  <- makeRateMatrix(2, 0.1)
  Q[1,2] <- 0.1
  Q[2,1] <- 0.2

  # this should succeed
  expect_silent(
    tp$setEtaConstantUnequal(Q)
  )

  # this should fail (negative)
  QQ <- makeRateMatrix(2, -0.1)
  expect_error(
    tp$setEtaConstantUnequal(QQ)
  )

  # this should fail (dimensions)
  QQ <- makeRateMatrix(3, 0.1)
  expect_error(
    tp$setEtaConstantUnequal( QQ )
  )

})

test_that("Setting functions: time-varying equal anagenetic matrix.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)

  # this should succeed
  expect_silent(
    tp$setEtaTimeVaryingEqual( 1, c(0.1, 0.2) )
  )

  # this should fail (time negative)
  expect_error(
    tp$setEtaTimeVaryingEqual( -1, c(0.1, 0.2) )
  )

  # this should fail (rate negative)
  expect_error(
    tp$setEtaTimeVaryingEqual( 1, c(-0.1, 0.2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setEtaTimeVaryingEqual(1, c(1,2,3) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setEtaTimeVaryingEqual( c(1,2), c(1,2) )
  )

})


test_that("Setting functions: time-varying equal anagenetic matrix.", {

  # create the dummy tensorphylo object
  tp <- new(TensorPhyloInstance, 2)
  Q  <- makeRateMatrix(2, 0.1)
  Q[1,2] <- 0.1
  Q[2,1] <- 0.2
  Qs <- c(Q, Q)

  # this should succeed
  expect_silent(
    tp$setEtaTimeVaryingUnequal( 1, Qs )
  )

  # this should fail (time negative)
  expect_error(
    tp$setEtaTimeVaryingUnequal( -1, Qs )
  )

  # this should fail (rate negative)
  QQ <- makeRateMatrix(2, -0.1)
  QQs <- c(QQ, QQ)
  expect_error(
    tp$setEtaTimeVaryingUnequal( 1, QQs )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setEtaTimeVaryingEqual(1, c(1,2,3) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setEtaTimeVaryingEqual( c(1,2), c(1,2) )
  )

  # this should fail (dimensions)
  expect_error(
    tp$setEtaTimeVaryingUnequal( c(1, 2), Qs  )
  )

})
