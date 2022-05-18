library(microbenchmark)
library(castor)
library(diversitree)
library(tensoRphylo)
library(TESS)
library(phangorn)
library(parallel)
library(stringr)

# parameter space
ntips   <- 2^(5:10)
nstates <- 2^(1:7)
nstates[nstates == 32] <- 31
reps    <- 1:100
calcs   <- 10

# generic settings
f      <- 1.0
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
  tree <- suppressWarnings(tess.sim.taxa(1, this_ntips, lambda = lambda, mu = mu, max = 10000)[[1]])
  # scale <- sum(tree$edge.length)
  # tree$edge.length <- tree$edge.length / scale
  # lambda <- lambda * scale
  # mu <- mu * scale
  # eta <- eta * scale

  Q <- matrix(eta / (this_nstates - 1), this_nstates, this_nstates)
  diag(Q) <- -eta
  data <- sim.history(tree, Q, message = FALSE)$states
  mode(data) <- "numeric"

  # make the tp likelihood
  data_tp <- matrix(0, length(data), this_nstates)
  rownames(data_tp) <- tree$tip.label
  colnames(data_tp) <- 1:this_nstates
  data_tp[cbind(1:length(data), data)] <- 1

  # simulate a random omega
  num_pairs <- this_nstates * this_nstates / 2
  omega <- vector("list", this_nstates)
  for(i in 1:this_nstates) {

    # make an empty matrix
    M <- matrix(0, this_nstates, this_nstates)

    u <- runif(1)
    M[i,i] <- u

    # get the number elements in the lower triangle
    k <- sum(lower.tri(M))

    this_element <- sample.int(k, size = 1)
    M[lower.tri(M)][this_element] <- 1 - u

    # store this slice
    omega[[i]] <- M

  }

  if ( this_nstates <= 31 ) {

    # get classe params
    pars <- diversitree:::starting.point.classe(tree, this_nstates)
    pars[grepl("lambda", names(pars))] <- 0
    pars[grepl("mu", names(pars))]     <- mu
    pars[grepl("q", names(pars))]      <- eta / (this_nstates - 1)

    # fill in the speciation rates
    for(i in 1:this_nstates) {

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
        if ( this_nstates > 9 ) {
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

    data[] <- NA
    classe_llf <- diversitree::make.classe(tree, data, k = this_nstates, sampling.f = rep(1, this_nstates), strict = FALSE)
    class_ll   <- classe_llf(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT)

  }

  tp <- makeTensorPhylo(tree, data_tp)
  tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setConditionalProbabilityType(conditionalProbability$TIME)
  tp$setSafeMode(FALSE)

  O <- makeCladogeneticEvents(this_nstates)
  for(i in 1:this_nstates) {

    # i is the ancestor state
    this_omega <- omega[[i]]

    # get the non-zero indices
    inds <- which(this_omega != 0, arr.ind = TRUE)

    # loop over inds
    for(j in 1:nrow(inds)) {

      # get the indices
      these_inds <- inds[j,]

      if ( these_inds[1] != these_inds[2] ) {
        # get the probability
        this_prob <- this_omega[these_inds[1], these_inds[2]]

        # assign value
        O[i,these_inds[1], these_inds[2]] <- this_prob
        # O[i,these_inds[2], these_inds[1]] <- this_prob * 0.5
      }

    }

  }

  # set parameters
  tp$setLambdaConstant(lambda)
  tp$setMuConstant(mu)
  tp$setEtaConstantEqual(eta)
  tp$setOmegaConstant(O)

  tensorphylo_llf <- function() tp$computeLogLikelihood()
  tensorphylo_ll  <- tensorphylo_llf()

  # do benchmarks
  # cat("  Benchmarking.\n")
  if ( this_nstates <= 31  ) {
    mp <- microbenchmark(
      tensorphylo = tensorphylo_llf(),
      diversitree = classe_llf(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT),
      times       = calcs,
      unit        = "milliseconds"
    )
  } else {
    mp <- microbenchmark(
      tensorphylo = tensorphylo_llf(),
      times       = calcs,
      unit        = "milliseconds"
    )
  }
  smp     <- summary(mp)
  means   <- smp$mean
  medians <- smp$median

  # don't forget to rm(tp) and gc()
  rm(tp); gc()

  # return
  res <- data.frame(ntips       = this_ntips,
                    nstates     = this_nstates,
                    rep         = this_rep,
                    likelihood  = tensorphylo_ll,
                    mean_time   = means[1],
                    median_time = medians[1],
                    method      = "tensorphylo")

  if ( this_nstates <= 31 ) {
    dt <- data.frame(ntips       = this_ntips,
                     nstates     = this_nstates,
                     rep         = this_rep,
                     likelihood  = class_ll,
                     mean_time   = means[2],
                     median_time = medians[2],
                     method      = "diversitree")
    res <- rbind(res, dt)
  } else {
    dt <- data.frame(ntips       = this_ntips,
                     nstates     = this_nstates,
                     rep         = this_rep,
                     likelihood  = NA,
                     mean_time   = NA,
                     median_time = NA,
                     method      = "diversitree")
    res <- rbind(res, dt)
  }

  return(res)

}, mc.preschedule = FALSE, mc.cores = 6))

# write the results
write.table(results, file = "results/classe_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
q()
