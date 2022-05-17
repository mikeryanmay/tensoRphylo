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
reps    <- 1:100
calcs   <- 10

# generic settings
f      <- 1
lambda <- 1.0 * f
mu     <- 0.5 * f
eta    <- 0.1 * f

# all combinations
all_combinations <- expand.grid(ntips = ntips, nstates = nstates, rep = reps)

# evaluate all combinations
results <- do.call(rbind, mclapply(1:nrow(all_combinations), function(i) {

  # get the settings
  row          <- all_combinations[i,]
  this_ntips   <- row[["ntips"]]
  this_nstates <- row[["nstates"]]
  this_rep     <- row[["rep"]]

  cat(this_ntips, " -- ", this_nstates, " -- ", this_rep, "\n", sep = "")

  # simulate the tree
  # cat("  Simulating tree.\n")
  tree <- suppressWarnings(tess.sim.taxa(1, this_ntips, lambda = lambda, mu = mu, max = 1000)[[1]])

  # make the parameters
  # this_eta <- eta * this_nstates
  Q <- matrix(eta / (this_nstates - 1), this_nstates, this_nstates)
  diag(Q) <- -eta

  # combine parameters for castor simulation
  params <- list(
    birth_rates = rep(lambda, this_nstates),
    death_rates = rep(mu, this_nstates),
    sampling_rates = rep(0, this_nstates),
    transition_matrix = Q
  )

  # simulate the character dataset
  # cat("  Simulating character dataset\n")
  sim  <- sim.history(tree, Q, message = FALSE)
  data <- sim$states
  mode(data) <- "numeric"

  data_dt <- matrix(1, nrow = length(tree$tip.label), ncol = this_nstates)
  rownames(data_dt) <- tree$tip.label
  colnames(data_dt) <- 1:this_nstates

  # compute the true likelihood
  treepar_ll  <- -as.numeric(TreePar::LikConstant(lambda, mu, 1, TreeSim::getx(tree), survival = 1))
  phydat      <- phyDat(t(t(data)), type = "USER", levels = 1:this_nstates)
  phangorn_ll <- pml(tree, phydat, rate = eta, bf = rep(1 / this_nstates, this_nstates))$logLik
  # true_ll     <- treepar_ll + phangorn_ll
  true_ll     <- treepar_ll

  # compute the likelihood with castor
  # cat("  Computing castor likelihood.\n")
  # castor_llf <- castor_musse_likelihood(tree, this_nstates, tip_pstates = data, root_prior = "flat", root_conditioning = "crown", verbose = FALSE)
  castor_llf <- castor_musse_likelihood(tree, this_nstates, tip_priors = data_dt, root_prior = "flat", root_conditioning = "crown", verbose = FALSE)
  castor_ll  <- castor_llf(params, 1)

  # compute the likelihood with diversitree
  # cat("  Computing diversitree likelihood.\n")
  q <- as.vector(Q)
  q <- q[q > 0]
  dt_pars <- c(rep(lambda, this_nstates), rep(mu, this_nstates), q)
  names(dt_pars) <- diversitree:::default.argnames.musse(this_nstates)

  data[] <- NA
  diversitree_llf <- make.musse(tree, data, this_nstates, strict = FALSE)
  diversitree_ll  <- diversitree_llf(dt_pars, condition.surv = TRUE, root=diversitree::ROOT.FLAT)

  # compute the likelihood with tensorphylo
  # cat("  Computing tensorphylo likelihood.\n")
  # data_tp <- matrix(0, nrow = length(tree$tip.label), ncol = this_nstates)
  # rownames(data_tp) <- tree$tip.label
  # colnames(data_tp) <- 1:this_nstates
  # data_tp[cbind(1:length(tree$tip.label), data)] <- 1.0

  # tpAuto <- makeTensorPhylo(tree, data_tp)
  tpAuto <- makeTensorPhylo(tree, data_dt)
  tpAuto$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tpAuto$setApplyTreeLikCorrection(FALSE)
  tpAuto$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  tpAuto$setLambdaConstant(lambda)
  tpAuto$setMuConstant(mu)
  tpAuto$setEtaConstantEqual(eta)
  tensorphyloAuto_llf <- function() tpAuto$computeLogLikelihood()
  tensorphyloAuto_ll  <- tensorphyloAuto_llf()

  # compare the likelihoods
  castor_error            <- abs(castor_ll             - true_ll)
  diversitree_error       <- abs(diversitree_ll        - true_ll)
  tensorphyloAuto_error   <- abs(tensorphyloAuto_ll    - true_ll)

  # do benchmarks
  # cat("  Benchmarking.\n")
  mp <- microbenchmark(
    castor            = castor_llf(params, 1),
    diversitree       = diversitree_llf(dt_pars, condition.surv = TRUE, root=diversitree::ROOT.FLAT),
    tensorphyloAuto   = tpAuto$computeLogLikelihood(),
    times             = calcs,
    unit              = "milliseconds"
  )
  smp     <- summary(mp)
  means   <- smp$mean
  medians <- smp$median

  # don't forget to rm(tp) and gc()
  rm(tp); gc()

  # return
  castor_df <- data.frame(ntips       = this_ntips,
                          nstates     = this_nstates,
                          rep         = this_rep,
                          true        = true_ll,
                          likelihood  = castor_ll,
                          error       = castor_error,
                          mean_time   = means[1],
                          median_time = medians[1],
                          method      = "castor")

  diversitree_df <- data.frame(ntips       = this_ntips,
                               nstates     = this_nstates,
                               rep         = this_rep,
                               true        = true_ll,
                               likelihood  = diversitree_ll,
                               error       = diversitree_error,
                               mean_time   = means[2],
                               median_time = medians[2],
                               method      = "diversitree")

  tensorphyloAuto_df <- data.frame(ntips       = this_ntips,
                                   nstates     = this_nstates,
                                   rep         = this_rep,
                                   true        = true_ll,
                                   likelihood  = tensorphyloAuto_ll,
                                   error       = tensorphyloAuto_error,
                                   mean_time   = means[3],
                                   median_time = medians[3],
                                   method      = "tensorphylo")

  res <- rbind(castor_df, diversitree_df, tensorphyloAuto_df)

  return(res)

}, mc.preschedule = FALSE, mc.cores = 6))

# write the results
write.table(results, file = "results/musse_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
q()
