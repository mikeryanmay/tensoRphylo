library(abind)
library(Matrix)

# read a tree
tree   <- read.nexus("test/tree.nex")
newick <- write.tree(tree)

# just a simple test
tp <- new(TensorPhylo, 2)
tp$setDebugMode(1)

tp$setLambdaConstant(-0.5)
tp$setLambdaConstant(-0.0)
tp$setLambdaConstant(0.0)
tp$setLambdaConstant(1.0)

tp$setLambdaStateVarying( c(1,2) )
tp$setLambdaStateVarying( c(-1,2) )

tp$setLambdaTimeVarying( 1 , c(1,2) )
tp$setLambdaTimeVarying( 1 , c(-1,2) )

tp$setLambdaTimeStateVarying( c(1.3), matrix(c(1,2,3,4), nrow = 2, byrow = TRUE) )

tp$setRootPrior( c(1,2) )
tp$setRootPrior( c(1,1) )
tp$setRootPrior( c(1,1) / 2 )
tp$setRootPrior( c(1,2) / 3 )

# tp$setRootPriorFlat()

tp$setRhoPresent(0.0)
tp$getRho()
tp$getRhoTimes()

tp$setRhoConstant(0.5, 0.5)

tp$setRhoConstant( c(0, 0.5), c(1, 0.5)  )

tp$setRhoStateVarying( c(0.5), matrix(c(1, 0.5), ncol = 2)  )

tp$setRhoStateVarying( c(0, 0.5), matrix(c(1, 0.5, 0.2, 0.7), ncol = 2)  )



tp$setXiConstant( 0, 0.5  )

tp$setXiStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )



tp$setUpsilonConstant( 0, 0.5  )

tp$setUpsilonStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )



tp$setEtaConstantEqual(0.5)

Q <- matrix(c(0,2,1,0), 2, 2)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)

tp$setEtaConstantUnequal(Q)

tp$setEtaTimeVaryingEqual( c(1), c(0.5, 1.0) )

M <- abind(Q, 2 * Q, along = 3)
tp$setEtaTimeVaryingUnequal( c(1), M  )



tp$setLambdaTimeVarying( c(1), c(0.5, 1.0) )

tp$setLambdaTimeStateVarying(
  c(1.3), matrix(c(1,2,3,4), nrow = 2, byrow = TRUE)
)


# lambda

tp$setLambdaConstant(1.2)

tp$setLambdaStateVarying( c(1,2) )






