library(microbenchmark)
library(rbenchmark)
library(castor)
library(diversitree)
library(tensoRphylo)
library(TESS)
library(phangorn)
library(parallel)

# parameter space
reps    <- 1:100

# all combinations
all_combinations <- expand.grid(rep = reps)

# evaluate all combinations
results <- do.call(rbind, lapply(1:nrow(all_combinations), function(i) {

  # get the settings
  this_rep <- all_combinations[i,]
  cat(this_rep, "\n", sep = "")

  repeat {

    # simulate parameter values
    lambda  <- rgamma(2, 2, 2)
    epsilon <- runif(1)
    mu      <- lambda * epsilon
    eta     <- rgamma(1, 2, rate = 4)

    # bisse parameters
    pars <- c(lambda, mu, eta, eta)

    # simulate the tree
    tree <- tree.bisse(pars, max.taxa = 500)
    if ( is.null(tree) == FALSE ) {
      break
    }

  }

  # get the data
  tree$node.label <- NULL
  data <- tree$tip.state
  data_tp <- matrix(0, nrow = length(tree$tip.label), ncol = 2)
  colnames(data_tp) <- c("0","1")
  rownames(data_tp) <- names(data)
  data_tp[cbind(1:nrow(data_tp), 1 + data)] <- 1.0

  # compute the likelihood with tensorphylo
  tp <- makeTensorPhylo(tree, data_tp)
  tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  tp$setLambdaStateDependent(lambda)
  tp$setMuStateDependent(mu)
  tp$setEtaConstantEqual(eta)
  tensorphylo_ll <- tp$computeLogLikelihood()

  # compute the likelihood wth diversitree
  bisse <- diversitree::make.bisse(tree, data)

  # compute the likelihood
  bisse_ll <- bisse(pars, root = 1, condition.surv = TRUE)

  # return
  res <- data.frame(rep           = this_rep,
                     x_likelihood = tensorphylo_ll,
                     y_likelihood = bisse_ll,
                     method       = "diversitree")

  return(res)

}))

# write the results
write.table(results, file = "results/bisse_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
# q()
