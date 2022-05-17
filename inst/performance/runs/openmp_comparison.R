library(microbenchmark)
library(rbenchmark)
library(castor)
library(diversitree)
library(tensoRphylo)
library(TESS)
library(phangorn)
library(parallel)
source("src/castor_likelihood.R")

# parameter space
ntips   <- 2^(5:10)
nstates <- 2^(1:5)
ncores  <- c(1,2,3,4)
reps    <- 1:100
calcs   <- 10

# generic settings
f      <- 1
lambda <- 1.0 * f
mu     <- 0.5 * f
eta    <- 0.1 * f

# all combinations
all_combinations <- expand.grid(ncores = ncores, nstates = nstates, ntips = ntips, rep = reps)

# evaluate all combinations
results <- do.call(rbind, lapply(1:nrow(all_combinations), function(i) {

  # get the settings
  row          <- all_combinations[i,]
  this_ntips   <- row[["ntips"]]
  this_nstates <- row[["nstates"]]
  this_ncores  <- row[["ncores"]]
  this_rep     <- row[["rep"]]

  cat(this_ntips, " -- ", this_nstates, " -- ", this_ncores, " -- ", this_rep, "\n", sep = "")

  # simulate the tree
  # cat("  Simulating tree.\n")
  tree <- suppressWarnings(tess.sim.taxa(1, this_ntips, lambda = lambda, mu = mu, max = 1000)[[1]])

  # make the parameters
  # this_eta <- eta * this_nstates
  Q <- matrix(eta / (this_nstates - 1), this_nstates, this_nstates)
  diag(Q) <- -eta

  # simulate the character dataset
  data_dt <- matrix(1, nrow = length(tree$tip.label), ncol = this_nstates)
  rownames(data_dt) <- tree$tip.label
  colnames(data_dt) <- 1:this_nstates

  # compute the likelihood with tensorphylo
  tp <- makeTensorPhylo(tree, data_dt)
  tp$setLikelihoodApproximator(approximatorVersion$PARALLEL_BRANCHWISE)
  tp$setNumberOfThreads(this_ncores)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setEtaConstantEqual(eta)
  mb <- microbenchmark(
    tp = tp$computeLogLikelihood(),
    times = calcs,
    unit  = "milliseconds"
  )
  smb     <- summary(mb)
  means   <- smb$mean
  medians <- smb$median

  # don't forget to rm(tp) and gc()
  rm(tp); gc()

  # return
  res <- data.frame(ntips       = this_ntips,
                    nstates     = this_nstates,
                    ncores      = this_ncores,
                    rep         = this_rep,
                    mean_time   = means,
                    median_time = medians,
                    method      = "tensorphylo")

  return(res)

}))

# write the results
write.table(results, file = "results/openmp_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
q()
