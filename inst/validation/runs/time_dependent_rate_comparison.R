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

  # simulate the tree
  # cat("  Simulating tree.\n")
  repeat {

    # simulate parameter values
    t0        <- 10
    lambda_1  <- rgamma(1, 1, 1)
    lambda_2  <- rgamma(1, 1, 1)
    t_shift   <- runif(1, 0, t0)
    epsilon   <- runif(1)
    mu        <- lambda_1 * epsilon
    lambda    <- function(t) {
      if ( t < t_shift ) {
        return(lambda_2)
      } else {
        return(lambda_2)
      }
    }

    # simulate tree
    tree <- suppressWarnings(tess.sim.age(1, lambda = lambda, mu = mu, age = t0, maxTaxa = 500)[[1]])
    if ( inherits(tree, "phylo") && length(tree$tip.label) > 3 ) {
      break
    }

  }

  # compute the likelihood with tensorphylo
  tp <- makeTensorPhylo(tree, num_states = 2)
  tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
  tp$setLambdaTimeDependent(t0 - t_shift, c(lambda_2, lambda_1))
  tp$setMuConstant(mu)
  tensorphylo_ll <- tp$computeLogLikelihood()

  # compute the likelihood with TESS
  tess_ll <- TESS::tess.likelihood.rateshift(branching.times(tree), c(lambda_1, lambda_2), mu, t_shift)

  # return
  tess <- data.frame(rep          = this_rep,
                     x_likelihood = tensorphylo_ll,
                     y_likelihood = tess_ll,
                     method = "TESS")

  return( tess )

}))

# write the results
write.table(results, file = "results/time_dependent_rate_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
# q()
