library(microbenchmark)
library(castor)
library(diversitree)
library(tensoRphylo)
library(TESS)
library(phangorn)
library(parallel)
library(stringr)
source("src/geohisse_likelihood.R")

# parameter space
ntips   <- 2^(5:10)
nhidden <- 0:9
reps    <- 1:100
calcs   <- 10

# all combinations
all_combinations <- expand.grid(ntips = ntips, nhidden = nhidden, rep = reps)

# generic settings
f      <- 1
lambda <- 1.0 * f
mu     <- 0.5 * f
eta    <- 0.1 * f

# evaluate all combinations
results <- do.call(rbind, mclapply(1:nrow(all_combinations), function(i) {

  # get the settings
  row          <- all_combinations[i,]
  this_ntips   <- row[["ntips"]]
  this_nhidden <- row[["nhidden"]]
  this_rep     <- row[["rep"]]

  cat(this_ntips, " -- ", this_nhidden, " -- ", this_rep, "\n", sep = "")

  # simulate the tree
  ## Generate a list with the parameters of the model:
  pars <- SimulateGeoHiSSE(hidden.traits = this_nhidden, return.GeoHiSSE_pars = TRUE)
  for(i in 1:(this_nhidden + 1)) {
    pars$model.pars[,i] <- c(lambda, lambda, lambda, mu, mu, eta, eta)
  }
  pars$q.01[upper.tri(pars$q.01)] <- eta
  pars$q.01[lower.tri(pars$q.01)] <- eta
  pars$q.0[upper.tri(pars$q.0)]   <- eta
  pars$q.0[lower.tri(pars$q.0)]   <- eta
  pars$q.1[upper.tri(pars$q.1)]   <- eta
  pars$q.1[lower.tri(pars$q.1)]   <- eta

  sink(tempfile())
  repeat {
    sim.geohisse <- SimulateGeoHiSSE(pars=pars, hidden.traits = this_nhidden, x0 = "01A", max.taxa = this_ntips)
    if (length(sim.geohisse$phy$tip.label) == this_ntips) {
      break
    }
  }
  sink(file = NULL)

  phy <- sim.geohisse$phy
  phy$node.label <- NULL
  sim.dat <- data.frame(taxon=sim.geohisse$data[,1], ranges=as.numeric(sim.geohisse$data[,2]))
  classe.pars <- sim.geohisse$classe.pars

  turnover <- rep(c(1,1,0), times = this_nhidden + 1)
  eps <- rep(c(1,1), times = this_nhidden + 1)
  trans.rate <- TransMatMakerGeoHiSSE(hidden.traits = this_nhidden)
  trans.rate.mod <- ParEqual(trans.rate, c(1,2))
  mod1 <- hisse::GeoHiSSE(phy = phy, data = sim.dat, f = c(1,1,1),
                          turnover = turnover, eps = eps,
                          hidden.states = this_nhidden > 0, trans.rate = trans.rate.mod,
                          turnover.upper = 100, trans.upper = 10,
                          starting.vals = c(lambda + mu, mu / lambda, eta))

  # make a tensorphylo classe model
  lambda_vec <- classe.pars[grepl("lambda", names(classe.pars))]
  mus        <- classe.pars[grepl("mu", names(classe.pars))]
  qs         <- classe.pars[grepl("q", names(classe.pars))]
  nstates    <- length(mus)

  dim <- 3 * (this_nhidden + 1)
  Q <- makeRateMatrix(dim)
  for(i in 1:dim) {
    for(j in 1:dim) {
      if (i != j) {
        if ( this_nhidden > 2 ) {
          # get the label
          Q[i,j] <- qs[paste0("q", paste0(str_pad(c(i,j), 2, pad = "0"), collapse = ""))]
        } else {
          Q[i,j] <- as.numeric(qs[paste0("q", i, j)])
        }
      }
    }
  }

  lambdas <- numeric(dim)
  W <- makeCladogeneticEvents(dim)
  for(i in 1:dim) {

    # get all the lambdas for this state
    if ( this_nhidden > 2  ) {
      these_lambdas <- lambda_vec[grepl(paste0("lambda", str_pad(i, 2, pad = "0")), names(lambda_vec))]
    } else {
      these_lambdas <- lambda_vec[grepl(paste0("lambda", i), names(lambda_vec))]
    }

    # compute the sum
    lambdas_sum <- sum(these_lambdas)
    lambdas[i]  <- lambdas_sum

    # normalize
    these_lambdas <- these_lambdas / lambdas_sum

    # set values in W
    c <- 1
    for(j in 1:dim) {
      for(k in j:dim) {
        if (j != k) {
          # get the label
          if ( this_nhidden > 2 ) {
            ancestor   <- str_pad(i, 2, pad = "0")
            these_inds <- str_pad(sort(c(j,k)), 2, pad = "0")
          } else {
            ancestor   <- str_pad(i, 1, pad = "0")
            these_inds <- str_pad(sort(c(j,k)), 1, pad = "0")
          }
          this_name <- paste0("lambda", paste0(c(ancestor, these_inds), collapse = ""))
          this_prob <- these_lambdas[this_name]
          if ( this_prob != 0.0 ) {
            W[i,j,k]  <- this_prob
          }
        }
      }
    }

  }

  # create the data
  dataA <- matrix(0, nrow = length(phy$tip.label), ncol = 3)
  dataA[cbind(1:nrow(dataA), sim.dat[,2] + 1)] <- 1.0
  colnames(dataA) <- paste0(c(0:2))
  rownames(dataA) <- phy$tip.label
  data <- do.call(cbind, lapply(1:(this_nhidden + 1), function(x) dataA))

  # make the tp object
  tp <- makeTensorPhylo(phy, data)
  tp$setApplyTreeLikCorrection(FALSE)
  tp$setLikelihoodApproximator(approximatorVersion$SEQUENTIAL_BRANCHWISE)
  tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

  tp$setLambdaStateVarying(lambdas)
  tp$setMuStateVarying(mus)
  tp$setEtaConstantUnequal(Q)
  tp$setOmegaConstant(W)
  # tp$computeLogLikelihood()

  mp <- microbenchmark(
    tensorphylo = tp$computeLogLikelihood(),
    hisse       = mod1(),
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
                    nhidden     = this_nhidden,
                    rep         = this_rep,
                    mean_time   = means[1],
                    median_time = medians[1],
                    method      = "tensorphylo")
  dt <- data.frame(ntips       = this_ntips,
                   nhidden     = this_nhidden,
                   rep         = this_rep,
                   mean_time   = means[2],
                   median_time = medians[2],
                   method      = "hisse")
  res <- rbind(res, dt)

  return(res)

}, mc.preschedule = FALSE, mc.cores = 6))

# write the results
write.table(results, file = "results/geohisse_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# quit
q()
