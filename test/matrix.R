library(expm)

# make a Q matrix
Q <- AnageneticMatrix(2, 0.1)

# try the assignment
tp <- new(TensorPhylo, 2)
tp$setDebugMode(tensoRphylo::debugMode$DBG_PRINT)

# a homogeneous model
tp$setEtaConstantUnequal(Q)

# a time-varying model
Qs <- c(Q, Q)
tp$setEtaTimeVaryingUnequal(1, Qs)


P <- TransitionMatrix(2)

P[1,2] <- 0.5
P[2,1] <- 0.1

Ps <- c(P, P)
