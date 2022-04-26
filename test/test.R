library(abind)
library(Matrix)

# read a tree
tree   <- read.nexus("test/tree.nex")
newick <- write.tree(tree)

# just a simple test
tp <- new(TensorPhylo, 2)

tp$setLambdaConstant(-0.5)
tp$setLambdaConstant(-0.0)
tp$setLambdaConstant(0.0)
tp$setLambdaConstant(1.0)

tp$setLambdaStateVarying( c(1,2) )
tp$setLambdaStateVarying( c(-1,2) )

tp$setLambdaTimeVarying( 1 , c(1,2) )
tp$setLambdaTimeVarying( 1 , c(-1,2) )
tp$getLambda()

tp$setLambdaTimeStateVarying( c(1.3), matrix(c(1,2,3,4), nrow = 2, byrow = TRUE) )
tp$getLambda()
tp$getLambdaTimes()

tp$getRootPrior()
tp$setRootPrior( c(1,2) )
tp$setRootPrior( c(1,1) )
tp$setRootPrior( c(1,1) / 2 )
tp$setRootPrior( c(1,2) / 3 )

# tp$setRootPriorFlat()

tp$getRho()
tp$getRhoTimes()

tp$setRhoPresent(0.0)
tp$getRho()
tp$getRhoTimes()

tp$setRhoConstant( 0.5, 0.5)
tp$getRho()
tp$getRhoTimes()

tp$setRhoConstant( c(0, 0.5), c(1, 0.5)  )
tp$getRho()
tp$getRhoTimes()

tp$setRhoStateVarying( c(0.5), matrix(c(1, 0.5), ncol = 2)  )
tp$getRho()
tp$getRhoTimes()

tp$setRhoStateVarying( c(0, 0.5), matrix(c(1, 0.5, 0.2, 0.7), ncol = 2)  )
tp$getRho()
tp$getRhoTimes()


tp$getXi()
tp$getXiTimes()

tp$setXiConstant( 0, 0.5  )
tp$getXi()
tp$getXiTimes()

tp$setXiStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )
tp$getXi()
tp$getXiTimes()



tp$getUpsilon()
tp$getUpsilonTimes()

tp$setUpsilonConstant( 0, 0.5  )
tp$getUpsilon()
tp$getUpsilonTimes()

tp$setUpsilonStateVarying( 0, matrix(c(0.5, 1.0), ncol = 2) )
tp$getUpsilon()
tp$getUpsilonTimes()



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






