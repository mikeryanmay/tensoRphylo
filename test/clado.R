library(microbenchmark)

# make an empty matrix
Omega <- CladogeneticEvents(num_states = 2)
Omega

# set some values
Omega[1,2,1] <- 0.25
Omega[1,1,2] <- 0.25
Omega

# try the assignment
tp <- new(TensorPhylo, 2)
tp$setDebugMode( tensoRphylo::debugMode$DBG_PRINT )



tp$setOmegaConstant(Omega)




Omegas <- c(Omega, Omega)

tp$setOmegaTimeVarying( 1, Omegas )

