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

tp$setRhoConstant(0.5, 0.5)

tp$setRhoConstant( c(0, 0.5), c(1, 0.5)  )

tp$setRhoStateVarying( c(0.5), matrix(c(1, 0.5), ncol = 2)  )

tp$setRhoStateVarying( c(0, 0.5), matrix(c(1, 0.5, 0.2, 0.7), ncol = 2)  )



tp$setXiConstant( 0, 0.5  )

tp$setXiStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )



tp$setUpsilonConstant( 0, 0.5  )

tp$setUpsilonStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )


tp$setGammaConstant( 0, 0.5  )

tp$setGammaStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )

P <- TransitionMatrix(2)
P[1,2] <- 0.5
P[2,1] <- 0.2

Ps <- c(P, P)

tp$setZeta( c(1), c(P) )



tp$setEtaConstantEqual(0.5)

Q <- AnageneticMatrix(2)
Q[1,2] <- 1
Q[2,1] <- 2

tp$setEtaConstantUnequal(Q)

tp$setEtaTimeVaryingEqual( c(1), c(0.5, 1.0) )


Qs <- c(Q, Q)

tp$setEtaTimeVaryingUnequal( c(1), Qs )




