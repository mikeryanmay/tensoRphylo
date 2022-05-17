library(tensoRphylo)
library(microbenchmark)

# read tree data
phy <- readRDS(system.file("testdata", "sampled_ancestor_tree.Rda", package = "tensoRphylo"))

# make parameters
lambda <- 0.1
mu     <- 0.03
phi    <- 0.05

# one core
tp <- makeTensorPhylo(phy, nstates = 100)

# set algorithm
tp$setLikelihoodApproximator(approximatorVersion$PARALLEL_BRANCHWISE)
tp$setNumberOfThreads(1)

# set parameters
tp$setLambdaConstant(lambda)
tp$setMuConstant(mu)
tp$setPhiConstant(phi)

# comparability: don't include the probability of the tree shape
tp$setApplyTreeLikCorrection(FALSE)

# compute
one_core <- microbenchmark(
  tp$computeLogLikelihood(),
  unit = "milliseconds"
)

# two cores
tp$setNumberOfThreads(2)
two_core <- microbenchmark(
  tp$computeLogLikelihood(),
  unit = "milliseconds"
)

# three cores
tp$setNumberOfThreads(3)
three_core <- microbenchmark(
  tp$computeLogLikelihood(),
  unit = "milliseconds"
)

# four cores
tp$setNumberOfThreads(4)
four_core <- microbenchmark(
  tp$computeLogLikelihood(),
  unit = "milliseconds"
)

# five cores
tp$setNumberOfThreads(5)
five_core <- microbenchmark(
  tp$computeLogLikelihood(),
  unit = "milliseconds"
)

tps <- vector("list", 4)
for(i in 1:4) {

  # make a tp instance
  tps[[i]] <- makeTensorPhylo(phy, nstates = 100)

  # set algorithm
  tps[[i]]$setLikelihoodApproximator(approximatorVersion$PARALLEL_BRANCHWISE)
  tps[[i]]$setNumberOfThreads(i)

  # set parameters
  tps[[i]]$setLambdaConstant(lambda)
  tps[[i]]$setMuConstant(mu)
  tps[[i]]$setPhiConstant(phi)

  # comparability: don't include the probability of the tree shape
  tps[[i]]$setApplyTreeLikCorrection(FALSE)

  cat(tps[[i]]$computeLogLikelihood(), "\n")

}

mb <- microbenchmark(
  tps[[1]]$computeLogLikelihood(),
  tps[[2]]$computeLogLikelihood(),
  tps[[3]]$computeLogLikelihood(),
  tps[[4]]$computeLogLikelihood(),
  unit = "milliseconds"
)

# # make a tp instance
# tp <- makeTensorPhylo(phy, nstates = 2)
#
# # set algorithm
# tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
# tp$setNumberOfThreads(1)
#
# # set parameters
# tp$setLambdaConstant(lambda)
# tp$setMuConstant(mu)
# tp$setPhiConstant(phi)
#
# # comparability: don't include the probability of the tree shape
# tp$setApplyTreeLikCorrection(FALSE)
#
# # compute
# one_core <- microbenchmark(
#   tp$computeLogLikelihood()
# )
#
# # more cores
# tp$setNumberOfThreads(4)
#
# four_core <- microbenchmark(
#   tp$computeLogLikelihood()
# )



