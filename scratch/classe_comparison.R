library(tensoRphylo)
library(diversitree)
library(TESS)
library(microbenchmark)
library(stringr)

# state space
num_states <- 31

# diversification rates
lambda <- 1.0
mu <- 0.5

# simulate a random omega
num_pairs   <- num_states * num_states / 2
omega <- vector("list", num_states)
for(i in 1:num_states) {

  # make an empty matrix
  M <- matrix(0, num_states, num_states)

  # get the number elements in the lower triangle
  k <- sum(lower.tri(M))

  # compute the number of non-zero elements
  num_nonzero <- ceiling(k * 0.1)

  # choose num_nonzero unique i,j at random
  these_elements <- sample.int(k, size = num_nonzero)

  # choose parameters
  probs <- rgamma(num_nonzero, 2, 1)
  probs <- probs / sum(probs)

  # assign parameters
  M[lower.tri(M)][these_elements] <- probs

  # store this slice
  omega[[i]] <- M

}

# eta
eta     <- 0.1
Q       <- matrix(eta / (num_states - 1), num_states, num_states)
diag(Q) <- -eta

# simulate some data
tree <- tess.sim.taxa(1, 100, 1000, lambda, mu)[[1]]
data <- sim.history(tree, Q, message = FALSE)$states
mode(data) <- "numeric"

# get classe params
pars <- diversitree:::starting.point.classe(tree, num_states)
pars[grepl("lambda", names(pars))] <- 0
pars[grepl("mu", names(pars))]     <- mu
pars[grepl("q", names(pars))]      <- eta / (num_states - 1)

# fill in the speciation rates
for(i in 1:num_states) {

  # i is the ancestor state
  this_omega <- omega[[i]]

  # get the non-zero indices
  inds <- which(this_omega != 0, arr.ind = TRUE)

  # loop over inds
  for(j in 1:nrow(inds)) {

    # get the indices
    these_inds <- inds[j,]

    # get the probability
    this_prob <- this_omega[these_inds[1], these_inds[2]]

    # build the name
    if ( num_states > 9 ) {
      ancestor   <- str_pad(i, 2, pad = "0")
      these_inds <- str_pad(sort(these_inds), 2, pad = "0")
    } else {
      ancestor   <- str_pad(i, 1, pad = "0")
      these_inds <- str_pad(sort(these_inds), 1, pad = "0")
    }
    this_name <- paste0("lambda", paste0(c(ancestor, these_inds), collapse = ""))

    # assign the parameter
    pars[this_name] <- lambda * this_prob

  }

}

# make the classe likelihood
model <- diversitree::make.classe(tree, data, k=num_states, sampling.f = rep(1, num_states), strict = FALSE)

# tensorphylo
data_tp <- matrix(0, length(data), num_states)
rownames(data_tp) <- tree$tip.label
colnames(data_tp) <- 1:num_states
data_tp[cbind(1:length(data), data)] <- 1

tp <- makeTensorPhylo(tree, data_tp)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$TIME)
tp$setSafeMode(FALSE)

O <- makeCladogeneticEvents(num_states)
for(i in 1:num_states) {

  # i is the ancestor state
  this_omega <- omega[[i]]

  # get the non-zero indices
  inds <- which(this_omega != 0, arr.ind = TRUE)

  # loop over inds
  for(j in 1:nrow(inds)) {

    # get the indices
    these_inds <- inds[j,]

    # get the probability
    this_prob <- this_omega[these_inds[1], these_inds[2]]

    # assign value
    O[i,these_inds[1], these_inds[2]] <- this_prob * 0.5
    O[i,these_inds[2], these_inds[1]] <- this_prob * 0.5

  }

}

# set parameters
tp$setLambdaConstant(lambda)
tp$setMuConstant(mu)
tp$setEtaConstantEqual(eta)
tp$setOmegaConstant(O)

# compare likelihoods
tp$computeLogLikelihood();model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT)

# compare runtimes
microbenchmark(
  dt = model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT),
  tp = tp$computeLogLikelihood(),
  times = 100
)

rm(tp)
gc(verbose = FALSE)














