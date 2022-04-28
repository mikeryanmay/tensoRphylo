library(expm)

# make a Q matrix
Q <- new(RateMatrix, 2, 0.1)
Q[1,2] <- 3.0

# try the assignment
tp <- new(TensorPhylo, 2)
tp$setDebugMode(tensoRphylo::debugMode$DBG_PRINT)

# a homogeneous model
tp$setEtaConstantUnequal(Q)

# a time-varying model
Qs <- c(Q, Q)
tp$setEtaTimeVaryingUnequal(1, Qs)

