library(microbenchmark)
library(rbenchmark)
library(castor)
library(diversitree)
library(tensoRphylo)
library(TESS)
library(phangorn)
library(parallel)

# parameter space
ntips   <- 2^(5:10)
nstates <- 2^(1:5)
algos   <- unlist(tensoRphylo::approximatorVersion[1:3])
reps    <- 1:100
calcs   <- 10

# generic settings
f      <- 1
lambda <- 1.0 * f
mu     <- 0.5 * f
eta    <- 0.1 * f

# all combinations
all_combinations <- expand.grid(ntips = ntips, nstates = nstates, algos = algos, rep = reps)

# evaluate all combinations
results <- do.call(rbind, mclapply(1:nrow(all_combinations), function(i) {

  # get the settings
  row            <- all_combinations[i,]
  this_ntips     <- row[["ntips"]]
  this_nstates   <- row[["nstates"]]
  this_algo      <- row[["algos"]]
  this_rep       <- row[["rep"]]
  this_algo_name <- names(tensoRphylo::approximatorVersion[this_algo + 1])

  cat(this_ntips, " -- ", this_nstates, " -- ", this_algo_name, " -- ", this_rep, "\n", sep = "")

  # simulate the tree
  # cat("  Simulating tree.\n")
  tree <- suppressWarnings(tess.sim.taxa(1, this_ntips, lambda = lambda, mu = mu, max = 1000)[[1]])

  # make the parameters
  Q <- matrix(eta / (this_nstates - 1), this_nstates, this_nstates)
  diag(Q) <- -eta

  # simulate the character dataset
  data_dt <- matrix(1, nrow = length(tree$tip.label), ncol = this_nstates)
  rownames(data_dt) <- tree$tip.label
  colnames(data_dt) <- 1:this_nstates

  # compute the true likelihood
  treepar_ll  <- -as.numeric(TreePar::LikConstant(lambda, mu, 1, TreeSim::getx(tree), survival = 1))
  true_ll     <- treepar_ll

  # make the tp object
  tp <- makeTensorPhylo(tree, data_dt)
  tp$setLikelihoodApproximator(this_algo)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
  tp$setIntegrationScheme(this_algo)

  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setEtaConstantEqual(eta)

  # compute the likelihood
  tp_ll <- tp$computeLogLikelihood()

  # do benchmarks
  mp <- microbenchmark(
    tensorphylo = tp$computeLogLikelihood(),
    times       = calcs,
    unit        = "milliseconds"
  )
  smp     <- summary(mp)
  means   <- smp$mean
  medians <- smp$median

  # don't forget to rm(tp) and gc()
  rm(tp); gc()

  # return
  res <- data.frame(ntips       = this_ntips,
                    nstates     = this_nstates,
                    this_algo   = this_algo,
                    error       = abs(true_ll - tp_ll),
                    rep         = this_rep,
                    mean_time   = means,
                    median_time = medians,
                    algo        = this_algo_name)

  return(res)

}, mc.cores = 6, mc.preschedule = FALSE))

# write the results
write.table(results, file = "results/approximator_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
q()
