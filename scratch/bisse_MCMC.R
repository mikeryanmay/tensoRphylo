library(tensoRphylo)
library(diversitree)

# tree
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)

# make likelihood
tree <- phy
tree$node.label <- NULL
data <- phy$tip.state

tp <- makeTensorPhylo(tree, data)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

lik <- function(pars) {

  # get parameters
  lambdas <- pars[1:2]
  mus     <- pars[3:4]
  Q       <- makeRateMatrix(2)
  Q[1,2]  <- pars[5]
  Q[2,1]  <- pars[6]

  # set parameters
  tp$setLambdaStateVarying(lambdas)
  tp$setMuStateVarying(mus)
  tp$setEtaConstantUnequal(Q)

  # return likelihood
  tp$computeLogLikelihood()

}

# lik <- make.bisse(tree, data)

# do MCMC?
tmp <- mcmc(lik, pars, nsteps=100, w=.1)
w <- diff(sapply(tmp[2:7], quantile, c(.05, .95)))
out <- mcmc(lik, pars, nsteps=1000, w=w)

profiles.plot(out["lambda0"], col.line="red")







