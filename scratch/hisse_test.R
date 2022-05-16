library(tensoRphylo)
library(diversitree)
library(hisse)
library(microbenchmark)

## Generate a list with the parameters of the model:
pars <- SimulateGeoHiSSE(hidden.traits = 1, return.GeoHiSSE_pars = TRUE)
pars

pars$model.pars[,1] <- c(0.1, 0.1, 0.1, 0.03, 0.03, 0.05, 0.05)
pars$model.pars[,2] <- c(0.2, 0.2, 0.2, 0.03, 0.03, 0.05, 0.05)
pars$q.01[1,2] <- pars$q.01[2,1] <- 0.005
pars$q.0[1,2] <- pars$q.0[2,1] <- 0.005
pars$q.1[1,2] <- pars$q.1[2,1] <- 0.005
pars

set.seed(42)
sim.geohisse <- SimulateGeoHiSSE(pars=pars, hidden.traits = 1, x0 = "01A", max.taxa = 500)

phy <- sim.geohisse$phy
phy$node.label <- NULL
sim.dat <- data.frame(taxon=sim.geohisse$data[,1], ranges=as.numeric(sim.geohisse$data[,2]))
classe.pars <- sim.geohisse$classe.pars

turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod1 <- hisse::GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
                 turnover = turnover, eps=eps,
                 hidden.states = FALSE, trans.rate = trans.rate.mod,
                 turnover.upper = 100, trans.upper = 10)
mod1()

# make a tensorphylo classe model
lambda_vec <- classe.pars[grepl("lambda", names(classe.pars))]
mus        <- classe.pars[grepl("mu", names(classe.pars))]
qs         <- classe.pars[grepl("q", names(classe.pars))]
nstates    <- length(mus)

Q <- makeRateMatrix(6)
for(i in 1:6) {
  for(j in 1:6) {
    if (i != j) {
      Q[i,j] <- as.numeric(qs[paste0("q", i, j)])
    }
  }
}

lambdas <- numeric(6)
W <- makeCladogeneticEvents(6)
for(i in 1:6) {

  # get all the lambdas for this state
  these_lambdas <- lambda_vec[grepl(paste0("lambda", i), names(lambda_vec))]

  # compute the sum
  lambdas_sum <- sum(these_lambdas)
  lambdas[i]  <- lambdas_sum

  # normalize
  these_lambdas <- these_lambdas / lambdas_sum

  # set values in W
  c <- 1
  for(j in 1:6) {
    for(k in j:6) {
      if (j != k) {
        c <- c + 1
        this_prob <- these_lambdas[c]
        if ( this_prob != 0.0 ) {
          W[i,j,k]  <- this_prob * 0.5
          W[i,k,j]  <- this_prob * 0.5
        }
      }
    }
  }

}

# create the data
dataA <- matrix(0, nrow = length(phy$tip.label), ncol = 3)
colnames(dataA) <- paste0(c(0:2),"A")
rownames(dataA) <- phy$tip.label
dataA[cbind(1:nrow(dataA), sim.dat[,2] + 1)] <- 1.0
dataB <- dataA
colnames(dataB) <- paste0(c(0:2),"B")
data  <- cbind(dataA, dataB)


# make the tp object
tp <- makeTensorPhylo(phy, data)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

tp$setLambdaStateVarying(lambdas)
tp$setMuStateVarying(mus)
tp$setEtaConstantUnequal(Q)
tp$setOmegaConstant(W)



microbenchmark(
  tp$computeLogLikelihood(),
  mod1()
)












