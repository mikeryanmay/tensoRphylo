library(abind)
library(Matrix)

# read a tree
tree   <- read.nexus("test/tree.nex")
newick <- write.tree(tree)

# just a simple test
tp <- new(TensorPhylo, 2)
# tp$report()

tp$getEta()

tp$setEtaConstantEqual(0.5)
tp$getEta()
tp$getEtaTimes()

Q <- matrix(c(0,2,1,0), 2, 2)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)

tp$setEtaConstantUnequal(Q)
tp$getEta()

tp$setEtaTimeVaryingEqual( c(1), c(0.5, 1.0) )
tp$getEta()

M <- abind(Q, 2 * Q, along = 3)
tp$setEtaTimeVaryingUnequal( c(1), M  )
tp$getEta()



tp$setLambdaTimeVarying( c(1), c(0.5, 1.0) )
tp$getLambda()
tp$getLambdaTimes()

tp$setLambdaTimeStateVarying(
  c(1.3), matrix(c(1,2,3,4), nrow = 2, byrow = TRUE)
)
tp$getLambda()
tp$getLambdaTimes()


tp$getMu()

# root frequency
tp$getRootPrior()


# lambda
tp$getLambda()

tp$setLambdaConstant(1.2)
tp$getLambda()

tp$setLambdaStateVarying( c(1,2) )
tp$getLambda()



# Q <- matrix(1, 2, 2)
# diag(Q) <- 0
# diag(Q) <- -rowSums(Q)
#
# # abind(list(Q, Q, Q), along = 3)
#
# M <- abind(Q, Q, along = 3)
# tp$eta <- M
# tp$eta
# tp$etaTimes
# tp$updateEtas()
#
# tp$eta <- array(Q, c(2,2,1))
# tp$eta
# tp$etaTimes
# tp$updateEtas()
#
#
# tp$lambda
# tp$lambdaTimes
# tp$updateLambdas()
#
# tp$mu
# tp$muTimes
# tp$updateMus()
#
# tp$phi
# tp$phiTimes
# tp$updatePhis()
#
# tp$delta
# tp$deltaTimes
# tp$updateDeltas()



# tp$lambda <- c(1.0, 2.0, 3.0)
# tp$lambda <- matrix(c(1.0, 2.0, 3.0))

# set the tree
# tp$setTree(newick)

# tp$setLambda( 1.0 )
# tp$setLambda( c(1.0, 2.0, 3.0) )
# tp$setLambda( c(1.0, 2.0, 3.0), c(1.0, 2.0, 3.0) )
# M <- matrix(0, 3, 3)
# tp$setLambda( c(1.0, 2.0, 3.0), M )
