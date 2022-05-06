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
k <- 2

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









tp <- new(TensorPhyloInstance, 2)

tp$setLambdaStateVarying(c(1,2))
tp$setMuConstant(0.1)
tp$setEtaConstantEqual(0.2)

O <- makeCladogeneticEvents(2)
tp$setOmegaConstant(O)

A <- tp$getQuasiStationaryFrequency(0)

# E <- eigen(A)
# E$vectors[,1] / sum(E$vectors[,1])

pars <- c(1,2,0.1,0.1,0.2,0.2)
diversitree:::stationary.freq.bisse(pars)










# read tree data
phy  <- tensoRphylo:::extant_tree
data <- tensoRphylo:::extant_data

# parameters
lambda <- 1.0
mu     <- 0.5
q      <- 0.1
pc     <- c(0.4, 0.3)
pa     <- c(0.7, 0.2)

# make a tp instance
tp <- makeTensorPhylo(phy, data)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$TIME)

# make the tp parameters
Q <- makeRateMatrix(2, q)
O <- makeCladogeneticEvents(2)
O[1,1,2] <- pc[1] * pa[1] * 0.5
O[1,2,1] <- pc[1] * pa[1] * 0.5
O[1,2,2] <- pc[1] * (1 - pa[1])
O[2,2,1] <- pc[2] * pa[2] * 0.5
O[2,1,2] <- pc[2] * pa[2] * 0.5
O[2,1,1] <- pc[2] * (1 - pa[2])

# set the parameters
tp$setLambdaConstant(lambda)
tp$setMuConstant(mu)
tp$setEtaConstantUnequal(Q)
tp$setOmegaConstant(O)

# compute the probability
tp$computeLogLikelihood()

# switch to classe params
pars <- diversitree::starting.point.classe(phy, 2)
pars["lambda111"] <- lambda * (1 - pc[1])
pars["lambda112"] <- lambda * pc[1] * pa[1]
pars["lambda122"] <- lambda * pc[1] * (1 - pa[1])
pars["lambda211"] <- lambda * pc[2] * (1 - pa[2])
pars["lambda212"] <- lambda * pc[2] * pa[2]
pars["lambda222"] <- lambda * (1 - pc[2])
pars["mu1"] <- mu
pars["mu2"] <- mu
pars["q12"] <- q
pars["q21"] <- q

data_vec  <- (data %*% c(0,1))[,1] + 1
model <- diversitree::make.classe(phy, data_vec, k=2, sampling.f=c(1,1), control = list(tol = 1e-16))
model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT)

# sprintf("%.10f", model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT))
# sprintf("%.10f", tp$computeLogLikelihood())


# E <- diversitree:::projection.matrix.classe(pars, 2)
# diversitree:::stationary.freq.classe.ev(pars, 2)
# diversitree:::projection.matrix.classe(pars, 2)
diversitree:::stationary.freq.classe.ev(pars, 2)
tp$getQuasiStationaryFrequency(0)




