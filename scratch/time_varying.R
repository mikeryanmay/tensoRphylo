library(tensoRphylo)
library(microbenchmark)
library(RevGadgets)
library(phytools)
library(TreePar)
library(TESS)

# read some stuff
tree <- read.nexus("test/tree.nex")

# make the instance
tp <- makeTensorPhylo(tree, nstates = 2)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType( conditionalProbability$ROOT_MRCA )

# make the model
lambda <- c(0.01, 0.02, 0.03)
mu     <- c(0.005, 0.01, 0.02)
t      <- c(25, 50)

# try tp
tp$setLambdaTimeVarying(t, lambda)
tp$setMuTimeVarying(t, mu)

tp$computeLogLikelihood()

# try TESS
times <- branching.times(tree)
tess.likelihood.rateshift(times, rev(lambda), rev(mu), rateChangeTimesLambda = max(times) - rev(t), rateChangeTimesMu = max(times) - rev(t))
