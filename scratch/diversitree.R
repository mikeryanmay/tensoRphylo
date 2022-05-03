library(diversitree)
library(microbenchmark)

# make the bisse tree
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
repeat {
  phy  <- tree.bisse(pars, max.t=30, x0=0)
  if ( is.null(phy) == FALSE ) {
    if ( length(unique(phy$tip.state)) > 1  ) {
      break
    }
  }
}

# make the bisse object
bisse <- make.bisse(phy, phy$tip.state)

# get the character data
dat <- phy$tip.state
data <- matrix(0, nrow = length(dat), ncol = 2)
rownames(data) <- names(dat)
colnames(data) <- 0:1
for(i in 1:length(dat)) {
  data[i,as.character(dat[i])] <- 1.0
}

# make tp
tp <- makeTensorPhylo(phy, data)
tp$setApplyTreeLikCorrection(FALSE)

tp$setLambdaStateVarying( pars[1:2] )
tp$setMuConstant( pars[3] )
tp$setEtaConstantEqual( pars[5] )


# compare the likelihoods
bisse(pars, root = ROOT.FLAT, condition.surv = FALSE)
tp$computeLogLikelihood()

mb <- microbenchmark(
  bisse(pars, root = ROOT.FLAT, condition.surv = FALSE),
  tp$computeLogLikelihood(),
  times = 1000L
)


# now try with musse
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32

repeat {
  phy <- tree.musse(pars, 30, x0=1)
  if ( is.null(phy) == FALSE ) {
    if ( length(unique(phy$tip.state)) > 2 ) {
      break
    }
  }
}

# make the musse object
musse <- make.musse(phy, phy$tip.state, 3)

# get the character data
dat  <- phy$tip.state
data <- matrix(0, nrow = length(dat), ncol = 3)
rownames(data) <- names(dat)
colnames(data) <- 1:3
for(i in 1:length(dat)) {
  data[i,as.character(dat[i])] <- 1.0
}

# make tp
tp <- makeTensorPhylo(phy, data)
tp$setApplyTreeLikCorrection(FALSE)

# make the parameters
tp$setLambdaStateVarying( c(.1,  .15,  .2) )
tp$setMuStateVarying(c(.03, .045, .06))

Q <- makeRateMatrix(3)
Q[1,2] <- 0.05
Q[1,3] <- 0
Q[2,1] <- 0.05
Q[2,3] <- 0.05
Q[3,1] <- 0
Q[3,2] <- 0.05

tp$setEtaConstantUnequal(Q)

# compare likelihoods
musse(pars, root = ROOT.FLAT, condition.surv = FALSE)
tp$computeLogLikelihood()

mb <- microbenchmark(
  musse(pars, root = ROOT.FLAT, condition.surv = FALSE),
  tp$computeLogLikelihood(),
  times = 100L
)

# now for bigger musse #
lambda <- rgamma(10, 2, 100)
mu     <- rgamma(10, 2, 200)
q      <- rep(0.01, 2 * choose(10,2) )

pars <- c(lambda, mu, q)

repeat {
  phy <- tree.musse(pars, 100, x0=1)
  if ( is.null(phy) == FALSE ) {
    if ( length(unique(phy$tip.state)) == 10 ) {
      break
    }
  }
}

# make the musse object
musse <- make.musse(phy, phy$tip.state, 10, control = list(tol = 1e-8))
musse(pars, root = ROOT.FLAT, condition.surv = FALSE)

# get the character data
dat  <- phy$tip.state
data <- matrix(0, nrow = length(dat), ncol = 10)
rownames(data) <- names(dat)
colnames(data) <- 1:10
for(i in 1:length(dat)) {
  data[i,as.character(dat[i])] <- 1.0
}


# make the tp object
tp <- makeTensorPhylo(phy, data)
tp$setApplyTreeLikCorrection(FALSE)

# make the parameters
tp$setLambdaStateVarying( lambda )
tp$setMuStateVarying( mu )

# Q <- makeRateMatrix(10, 9 * 0.01)
# tp$setEtaConstantUnequal( Q )
# tp$computeLogLikelihood()

tp$setEtaConstantEqual( 9 * 0.01 )
tp$computeLogLikelihood()


mb <- microbenchmark(
  musse(pars, root = ROOT.FLAT, condition.surv = FALSE),
  tp$computeLogLikelihood(),
  times = 100L
)

x1 <- musse(pars, root = ROOT.FLAT, condition.surv = FALSE)
x2 <- tp$computeLogLikelihood()

sprintf("%.10f", x1)
sprintf("%.10f", x2)





