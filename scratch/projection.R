library(tensoRphylo)
library(diversitree)
library(expm)

# parameters
lambda <- 1
mu     <- 0.5
q      <- c(0.1,0.2)
pc     <- c(0.4, 0.3)
pa     <- c(0.7, 0.2)

# switch to classe params
phy <- tensoRphylo:::extant_tree

pars <- diversitree::starting.point.classe(phy, 2)
pars["lambda111"] <- lambda * (1 - pc[1])
pars["lambda112"] <- lambda * pc[1] * pa[1]
pars["lambda122"] <- lambda * pc[1] * (1 - pa[1])
pars["lambda211"] <- lambda * pc[2] * (1 - pa[2])
pars["lambda212"] <- lambda * pc[2] * pa[2]
pars["lambda222"] <- lambda * (1 - pc[2])
pars["mu1"] <- mu
pars["mu2"] <- mu
pars["q12"] <- q[1]
pars["q21"] <- q[2]

# make the projection matrix function
projection.matrix.classe <- function (pars, k) {
  A <- matrix(0, nrow = k, ncol = k)
  nsum <- k * (k + 1)/2
  kseq <- seq_len(k)
  pars.lam <- pars[seq(1, nsum * k)]
  pars.mu <- pars[seq(nsum * k + 1, (nsum + 1) * k)]
  pars.q <- pars[seq((nsum + 1) * k + 1, length(pars))]
  idx.lam <- cbind(rep(kseq, each = nsum), rep(rep(kseq, times = seq(k, 1, -1)), k), unlist(lapply(kseq, function(i) i:k)))
  idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])),
                 rep(kseq, each = k - 1))
  for (n in seq_len(nsum * k)) {
    r <- idx.lam[n, ]
    A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
    A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
  }
  A[idx.q] <- A[idx.q] + pars.q
  diag(A) <- 0
  diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) sum(pars.lam[seq((i - 1) * nsum + 1, i * nsum)]) - pars.mu[i]))
  A
}

A <- projection.matrix.classe(pars, 2)

diversitree:::stationary.freq.classe.ev(pars, k)
diversitree:::stationary.freq.classe(pars, k)
















