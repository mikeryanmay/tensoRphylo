library(tensoRphylo)
library(microbenchmark)
library(RevGadgets)
library(phytools)
library(TreePar)

# read some stuff
tree <- read.nexus("test/tree.nex")

# create the tp object
tp <- makeTensorPhylo(tree, nstates = 2)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType( conditionalProbability$ROOT_MRCA )

# make a likelihood function
obj <- function(pars) {

  lambda <- exp(pars[1])
  mu     <- exp(pars[2])

  # set the values
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)

  ll <- tp$computeLogLikelihood()

  cat(ll, lambda, mu, "\n", sep = "\t")

  return(-ll)

}

# initial values
init <- log(c(0.2, 0.05))

obj(init)

# fit the model
fit <- optim(init, obj, control = list(maxit = 1000))
fit$par <- exp(fit$par)

lambda <- fit$par[1]
mu     <- fit$par[2]

tp$setLambdaConstant(lambda)
tp$setMuConstant(mu)

tp$computeLogLikelihood()
LikConstant(lambda, mu, 1, getx(tree), survival = 0)



