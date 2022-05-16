library(microbenchmark)
library(castor)
library(diversitree)
library(tensoRphylo)
library(phangorn)
source("src/castor_likelihood.R")

# make the parameters
nstates <- 3
r <- 0.1
Q <- matrix(r / (nstates - 1), nstates, nstates)
diag(Q) <- -r
lambda <- 1
mu <- 0.5

params <- list(
  birth_rates = rep(lambda, nstates),
  death_rates = rep(mu, nstates),
  transition_matrix = Q
)

for(i in 1:10) {

  cat(i, "\n", sep="")

  # simulate a tree
  repeat {

    sim  <- castor::simulate_musse(nstates,
                                   parameters = params,
                                   max_time = 10,
                                   include_labels = TRUE,
                                   tip_basename = "t_")
    tree <- sim$tree
    data <- sim$tip_states

    # make sure every state is observed
    if (length(unique(data)) == nstates) {
      break
    }

  }

  # format the data for tp
  data_tp <- matrix(0, nrow = length(tree$tip.label), ncol = nstates)
  rownames(data_tp) <- tree$tip.label
  colnames(data_tp) <- 1:nstates
  data_tp[cbind(1:length(tree$tip.label), data)] <- 1.0

  # compute probability with tp
  tp <- makeTensorPhylo(tree, data_tp)
  tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setEtaConstantEqual(r)
  tp$computeLogLikelihood()

  # compute the probability with castor
  castor_llf <- castor_musse_likelihood(tree, nstates, tip_pstates = data, root_prior = "flat", root_conditioning = "crown")
  castor_llf(params, 1)

  # compute the probability with diversitree
  q <- as.vector(Q)
  q <- q[q > 0]
  dt_pars <- c(rep(lambda, nstates), rep(mu, nstates), q)
  names(dt_pars) <- diversitree:::default.argnames.musse(nstates)

  musse_llf <- make.musse(tree, data, nstates)
  musse_llf(dt_pars, condition.surv = TRUE, root=diversitree::ROOT.FLAT)

  mp <- microbenchmark(
    tensoRphylo = tp$computeLogLikelihood(),
    diversitree = musse_llf(dt_pars, condition.surv = FALSE, root=diversitree::ROOT.FLAT),
    castor = castor_llf(params, 1),
    times = 100
  )

}

q()
