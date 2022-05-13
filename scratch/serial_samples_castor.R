library(ape)
library(tensoRphylo)
library(castor)
source("inst/performance/src/castor_likelihood.R")

# read stuff
# tree <- read.nexus("inst/testdata/sampled_ancestor_tree.nex")
# data <- readNexusData("inst/testdata/sampled_ancestor_data.nex")

tree <- readRDS("inst/testdata/serial_sampled_tree.Rda")
data <- readRDS("inst/testdata/serial_sampled_data.Rda")

# remove the root edge
tree$root.edge <- NULL

# turn data into a vector
data_vec  <- 1 + (data %*% c(0,1))[,1]

# make the likelihood functions
lambda <- 0.1
mu     <- 0.03
delta  <- 0.05
eta    <- 0.1

Q <- matrix(eta / (2 - 1), 2, 2)
diag(Q) <- -eta

# combine parameters for castor simulation
params <- list(
  birth_rates = rep(lambda, 2),
  death_rates = rep(mu, 2),
  sampling_rates = rep(delta, 2),
  transition_matrix = Q
)

castor_llf <- castor_musse_likelihood(tree, 2, sampling_rates = delta, tip_pstates = data_vec, root_prior = "flat", root_conditioning = "crown", verbose = FALSE)
castor_llf(params, 1)

# tensorphylo

# make a tp instance
tp <- makeTensorPhylo(tree, data)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$ROOT_SAMPLING_AND_MRCA)

# set parameters
tp$setLambdaConstant(lambda)
tp$setMuConstant(mu)
tp$setDeltaConstant(delta)
tp$setEtaConstantEqual(eta)

# compute probability
tp$computeLogLikelihood()




