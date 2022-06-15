library(TreePar)
library(Rcpp)
library(dispRity)
sourceCpp("src/bdss_functions.cpp")

tree_logprob <- function(branching_times, sampling_times, root = 0, survival = T, lambda, mu, psi, rho = 0) {

  nbirthevents <- length(branching_times)
  nsampletips  <- length(sampling_times)
  if ( root == 0 ) {
    nextanttips <- nbirthevents - nsampletips
  } else {
    nextanttips <- nbirthevents + 1 - nsampletips
  }

  lik <- 0
  # rho^nextanttips
  if (nextanttips > 0) {
    if (rho == 0) {
      stop("There should not exist any extant tips when the sampling probability at the present is zero.")
    } else {
      lik <- lik + nextanttips * log(rho)
    }
  }

  # lambda^nbirthevents / (q(x_0)*q(x_1)*...*q(x_nbirthevents))
  logqxs <- numeric(nbirthevents)
  for (i in 1:nbirthevents) {
    logqxs[i] <- qfuncC(branching_times[i], lambda, mu, psi, rho, logprob = T)
  }
  lik <- lik + (nbirthevents - (root + 1)) * log(lambda) - sum(logqxs)

  # psi^nsampletips * (q(y_0)*q(y_1)*...*q(y_nsampletips))
  # here we assume that a sampling event induce an immediate death (even if the virus lineage is still alive it can never transmit any more)
  # i.e., we assume r = 1
  if (nsampletips > 0) {
    logqys <- numeric(nsampletips)
    for (i in 1:nsampletips) {
      logqys[i] <- qfuncC(sampling_times[i], lambda, mu, psi, rho, logprob = T)
    }
    lik <- lik + nsampletips * log(psi) + sum(logqys)
  }

  if (survival) { # condition on survival or not
    lik <- lik - (root + 1) * log(1 - p0sersampC(max(branching_times), lambda, mu, psi, rho))
  }

  return(lik)
}
