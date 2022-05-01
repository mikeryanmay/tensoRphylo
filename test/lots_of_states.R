library(tensoRphylo)
library(microbenchmark)
library(phangorn)

# numeric stuff
NUM_CATS = 5

# read some stuff
tree <- read.nexus("test/primates_tree.nex")
tp   <- tensoRphylo(tree, nstates = NUM_CATS * NUM_CATS)

# set the conditional probability type
# tp$setConditionalProbabilityType( tensoRphylo::conditionalProbability$ROOT_MRCA )

# set the sampling fraction
tp$setRhoPresent( 233 / 367 )

# make an objective function
obj <- function(pars) {

  # get discretized speciation rates
  lambda_mean  <- pars[1]
  lambda_alpha <- pars[2]
  lambdas      <- lambda_mean * discrete.gamma(lambda_alpha, NUM_CATS)
  if ( any(lambdas < 0.0) ) {
    return(-Inf)
  }

  # get discretized extinction rates
  mu_mean  <- pars[3]
  mu_alpha <- pars[4]
  mus      <- mu_mean * discrete.gamma(mu_alpha, NUM_CATS)
  if ( any(mus < 0.0) ) {
    return(-Inf)
  }

  # if ( mu_mean > lambda_mean  ) {
  #   return(-Inf)
  # }

  # get transition rates
  eta <- pars[5]
  if ( any(eta < 0.0) ) {
    return(-Inf)
  }

  # combine the diversification rates
  lambda_comb <- rep(lambdas, each = NUM_CATS)
  mu_comb     <- rep(mus, times = NUM_CATS)

  # set the values
  tp$setLambdaStateVarying(lambda_comb)
  tp$setMuStateVarying(mu_comb)
  tp$setEtaConstantEqual(eta)

  ll <- tp$computeLogLikelihood()

  cat(ll, lambda_mean, lambda_alpha, mu_mean, mu_alpha, eta, "\n", sep = "\t")

  return(-ll)

}

# initial values
init <- c(0.2, 2, 0.05, 2, 0.001)

# fit the model
fit <- optim(init, obj, control = list(maxit = 5000), method = "Nelder-Mead")

# plot the rates
curve(dgamma(x, fit$par[2], fit$par[2] / fit$par[1]), from = 0, to = 4 * fit$par[1], n = 1001, col = "blue")
curve(dgamma(x, fit$par[4], fit$par[4] / fit$par[3]), n = 1001, add = TRUE, col = "red")

# obj( c(0.15, 1.6, 0.7, 2.2, 0.007) )

# stochastic maps
branch_rates     <- tp$drawBranchRates(100)
avg_branch_rates <- Reduce("+", branch_rates) / length(branch_rates)

# make the color ramp
limits <- range(pretty(range(avg_branch_rates[,1])))
bins   <- seq(limits[1], limits[2], length.out=101)
cols   <- colorRampPalette(c("blue","green"))(100)
branch_cols <- cols[findInterval(avg_branch_rates[,1], bins)]

pdf("test/BDS.pdf", height = 5, width = 5)
plot.phylo(tree, edge.col = branch_cols, show.tip.label = FALSE, edge.width = 2, no.margin = TRUE)
dev.off()










