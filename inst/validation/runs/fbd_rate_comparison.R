library(diversitree)
library(tensoRphylo)
library(TESS)
library(TreeSim)
library(phangorn)
library(parallel)
library(FossilSim)
source("src/bdss_prob.R")

# parameter space
reps    <- 1:100

# all combinations
all_combinations <- expand.grid(rep = reps)

# evaluate all combinations
results <- do.call(rbind, lapply(1:nrow(all_combinations), function(i) {

  # get the settings
  this_rep <- all_combinations[i,]
  cat(this_rep, "\n", sep = "")

  # simulate parameter values
  lambda  <- rgamma(1, 1, 1)
  epsilon <- runif(1)
  mu      <- lambda * epsilon
  phi     <- lambda * runif(1)

  # simulate the tree
  # cat("  Simulating tree.\n")
  repeat {
    tree <- sim.fbd.age(log(50) / (lambda - mu), 1, lambda, mu, phi)[[1]]
    SAtree(tree)
    beast.fbd.format(tree)
    if ( inherits(tree, "phylo") ) {
      break
    }
  }

  # simulate fossils

  # compute the likelihood with tensorphylo
  tp <- makeTensorPhylo(tree, num_states = 2)
  tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tensorphylo_ll <- tp$computeLogLikelihood()

  # compute the likelihood with TreePar
  treepar_ll <- -as.numeric(TreePar::LikConstant(lambda, mu, 1, TreeSim::getx(tree), survival = 1))

  # compute the likelihood wth TESS
  tess_ll <- TESS::tess.likelihood(branching.times(tree), lambda = lambda, mu = mu)

  # return
  tess <- data.frame(rep          = this_rep,
                     x_likelihood = tensorphylo_ll,
                     y_likelihood = tess_ll,
                     method = "TESS")

  treepar <- data.frame(rep          = this_rep,
                        x_likelihood = tensorphylo_ll,
                        y_likelihood = treepar_ll,
                        method = "TreePar")

  return( rbind(tess, treepar) )

}))

# write the results
write.table(results, file = "results/constant_rate_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
# q()
