library(tensoRphylo)
library(TESS)

# make parameters
lambda <- 0.1
mu     <- 0.03
eta    <- 0.05
tree   <- tess.sim.taxa(1, 50, lambda, mu, max = 1000)[[1]]

tps <- vector("list", 100)
for(i in 1:100) {
  tps[[i]] <- makeTensorPhylo(tree, nstates = 2 + 1)
  tps[[i]]$setApplyTreeLikCorrection(FALSE)
  tps[[i]]$setConditionalProbabilityType(conditionalProbability$TIME)
  tps[[i]]$setLambdaConstant(lambda)
  tps[[i]]$setMuConstant(mu)
  tps[[i]]$setEtaConstantEqual(eta)
  tps[[i]]$computeLogLikelihood()
  cat(i, "\n")
}

# tps[[2]]$setLambdaConstant(0.2)
# tps[[2]]$computeLogLikelihood()
# tps[[1]]$computeLogLikelihood()

sapply(tps, function(x) x$computeLogLikelihood())

# # make a tp instance
# tpAuto <- makeTensorPhylo(tree, data)
# tpAuto$setApplyTreeLikCorrection(FALSE)
# tpAuto$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
# tpAuto$setLambdaConstant(lambda)
# tpAuto$setMuConstant(mu)
# tpAuto$setEtaConstantEqual(eta)
# tpAuto$computeLogLikelihood()
#
# tpBranch <- makeTensorPhylo(tree, data)
# tpBranch$setApplyTreeLikCorrection(FALSE)
# tpBranch$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
# tpBranch$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
# tpBranch$setLambdaConstant(lambda)
# tpBranch$setMuConstant(mu)
# tpBranch$setEtaConstantEqual(eta)
# tpBranch$computeLogLikelihood()
#
# tpSync <- makeTensorPhylo(tree, data)
# tpSync$setApplyTreeLikCorrection(FALSE)
# tpSync$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_OPTIMIZED)
# tpSync$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
# tpSync$setLambdaConstant(lambda)
# tpSync$setMuConstant(mu)
# tpSync$setEtaConstantEqual(eta)
# tpSync$computeLogLikelihood()
#
# tpAuto$computeLogLikelihood()
# tpSync$computeLogLikelihood()
# # tpBranch$computeLogLikelihood()
