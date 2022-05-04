library(bdms)
library(geiger)
library(tensoRphylo)
library(microbenchmark)
source("scratch/fbd_equations.R")
source("scratch/bdss_prob.R")

# read some data
tree <- read.nexus("scratch/tree_bdms.nex")
seq  <- read.nexus.data("scratch/seq_bdms.nex")

# drop the first site
seq <- do.call(rbind, seq)
# seq <- matrix(seq, ncol = 1)

# parameters
lambda <- 0.2
delta  <- 0.0
phi    <- 0.1
gamma  <- 0.05

# make the yule model
model <- rep("-", length = 10)
# model[1] <- "A"
yule_func <- YuleLikelihood2(tree, seq, model = model, lambda, gamma, delta, 0, phi)

# make the TP model
tp <- makeTensorPhylo(tree, nstates = 4)
tp$setApplyTreeLikCorrection(FALSE)

# tp parameters
lambdas <- c(lambda + delta, lambda, lambda, lambda)
deltas  <- phi
eta     <- gamma

tp$setLambdaConstant(lambda)
tp$setDeltaConstant(phi)
tp$setEtaConstantEqual(eta)

# J's function
ages         <- tree.age(tree, digits = 6)
branch_times <- ages[-c(1:length(tree$tip.label)),1]
branch_times <- c(branch_times, max(branch_times) + tree$root.edge)
sample_times <- ages$ages[ages$ages > 0]
sample_times <- sample_times[sample_times %in% branch_times == FALSE]

# loop over lambdas
lambda_x <- seq(0.1, 1, 0.01)
lambda_probs_yule <- numeric(length(lambda_x))
lambda_probs_tp   <- numeric(length(lambda_x))
lambda_probs_J    <- numeric(length(lambda_x))
for(i in 1:length(lambda_x)) {

  # get the value of lambda
  this_lambda <- lambda_x[i]

  # set values
  yule_func$setLambda0(this_lambda)
  tp$setLambdaConstant(this_lambda)

  # compute likelihoods
  yule_func$computeLikelihood()

  lambda_probs_yule[i] <- yule_func$selected_site_log_likelihood
  lambda_probs_tp[i]   <- tp$computeLogLikelihood()
  lambda_probs_J[i]    <- tree_logprob(branch_times, sample_times, root = 0, survival = TRUE, this_lambda, 0, phi, rho = 1)

}

plot(lambda_probs_yule - lambda_probs_tp)
plot(lambda_probs_yule - lambda_probs_J)
plot(lambda_probs_J    - lambda_probs_tp)

pdf("~/Downloads/likelihood_surface_lambda.pdf", height = 5, width = 10)
par(mar=c(4,4,0,0)+0.1)
plot(lambda_x, lambda_probs_yule, type = "p", lwd = 1, col = "black", xlab = "lambda", ylab = "likelihood", cex = 1)
lines(lambda_x, lambda_probs_J,   type = "p", pch = 3, lwd = 1, cex = 1)
lines(lambda_x, lambda_probs_tp,  type = "p", pch = 4, lwd = 1, lty = 1, cex = 1)
legend("topright", c("YuleSelection","J's function", "tensorphylo"), pch = c(1,3,4), bty = "n")
dev.off()

# loop over sampling rates
yule_func$setLambda0(lambda)
tp$setLambdaConstant(lambda)

phi_x <- seq(0.1, 1, 0.01)
phi_probs_yule <- numeric(length(phi_x))
phi_probs_tp   <- numeric(length(phi_x))
phi_probs_J    <- numeric(length(phi_x))

for(i in 1:length(phi_x)) {

  # get the value of phi
  this_phi <- phi_x[i]

  # set values
  yule_func$setPhi(this_phi)
  tp$setDeltaConstant(this_phi)

  # compute likelihoods
  yule_func$computeLikelihood()

  phi_probs_yule[i] <- yule_func$selected_site_log_likelihood
  phi_probs_tp[i]   <- tp$computeLogLikelihood()
  phi_probs_J[i]    <- tree_logprob(branch_times, sample_times, root = 0, survival = TRUE, lambda, 0, this_phi, rho = 1)

}

pdf("~/Downloads/likelihood_surface_phi.pdf", height = 5, width = 10)
par(mar=c(4,4,0,0)+0.1)
plot(phi_x, phi_probs_yule, type = "p", lwd = 1, col = "black", xlab = "phi", ylab = "likelihood", cex = 1)
lines(phi_x, phi_probs_J,   type = "p", pch = 3, lwd = 1, cex = 1)
lines(phi_x, phi_probs_tp,  type = "p", pch = 4, lwd = 1, lty = 1, cex = 1)
legend("topright", c("YuleSelection","J's function", "tensorphylo"), pch = c(1,3,4), bty = "n")
dev.off()


# now with selection


# yule_func$computeSelectedLikelihood()
# yule_func$selected_site_log_likelihood
# yule_func$computeNeutralLikelihood()
# yule_func$neutral_site_log_likelihoods

# now with tensorphylo
dat <- t(t(seq[,1]))
data <- matrix(0, nrow(dat), 4)
colnames(data) <- c("a","c","g","t")
rownames(data) <- rownames(dat)
for(i in 1:nrow(dat)) {
  data[i, dat[i,1]] <- 1.0
}

lambdas <- c(lambda + delta, lambda, lambda, lambda)
deltas  <- phi
eta     <- gamma # ?

# tp <- makeTensorPhylo(tree, data)
tp <- makeTensorPhylo(tree, nstates = 4)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$TIME)

tp$setLambdaStateVarying(lambdas)
tp$setDeltaConstant(deltas)
tp$setEtaConstantEqual(eta)

tp$computeLogLikelihood()
yule_func$selected_site_log_likelihood





n_tips    <- length(tree$tip.label)
n_extant  <- length(drop.extinct(tree)$tip.label)
n_extinct <- n_tips - n_extant

branch_prob <- dpois(0, (lambda + delta) * sum(tree$edge.length), log = TRUE)
tip_prob    <- n_extinct * log(phi)
node_part   <- tree$Nnode * log(lambda)

branch_prob + tip_prob + node_part

tree_logprob(branch_times, sample_times, root = 0, survival = TRUE, lambda, 0, phi, rho = 1)




# now with selection

# parameters
speciation_rate <- 0.2
selection       <- 0.1
sampling_rate   <- 0.1
mutation_rate   <- 0.05

# make the yule model
model[1] <- "A"
yule_func <- YuleLikelihood2(tree, t(t(seq[,1])), model = "A", speciation_rate, mutation_rate, selection, 0, sampling_rate)
yule_func$computeLikelihood()

# make the tp model
dat <- t(t(seq[,1]))
data <- matrix(0, nrow(dat), 4)
colnames(data) <- c("a","c","g","t")
rownames(data) <- rownames(dat)
for(i in 1:nrow(dat)) {
  data[i, dat[i,1]] <- 1.0
}

tp <- makeTensorPhylo(tree, data)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$TIME)

tp$setLambdaStateVarying( c(speciation_rate + selection, speciation_rate, speciation_rate, speciation_rate) )
tp$setDeltaConstant( sampling_rate )
tp$setEtaConstantEqual( mutation_rate )

# try different deltas
deltas <- seq(0.01, 0.1, 0.001)
ys_probs <- numeric(length(deltas))
tp_probs <- numeric(length(deltas))
for(i in 1:length(deltas)) {

  # get the delta
  this_delta <- deltas[i]

  # set the values
  yule_func$setDelta(this_delta)
  tp$setLambdaStateVarying( c(speciation_rate + this_delta, speciation_rate, speciation_rate, speciation_rate) )

  # compute the likelihoods
  ys_probs[i] <- yule_func$computeLikelihood()
  tp_probs[i] <- tp$computeLogLikelihood()

}

pdf("~/Downloads/likelihood_surface_delta.pdf", height = 5, width = 10)
par(mar=c(4,4,0,0)+0.1)
plot(phi_x, ys_probs, type = "p", lwd = 1, col = "black", xlab = "delta", ylab = "likelihood", cex = 1)
lines(phi_x, tp_probs,  type = "p", pch = 4, lwd = 1, lty = 1, cex = 1)
legend("topleft", c("YuleSelection", "tensorphylo"), pch = c(1,3,4), bty = "n")
dev.off()


# yule_func$setDelta(sampling_rate)
# tp$setDeltaConstant(sampling_rate)
#
# ys_foo <- function() {
#   yule_func$sel_lik_dirty <- TRUE
#   yule_func$computeLikelihood()
# }
#
# microbenchmark(
#   ys_foo(),
#   tp$computeLogLikelihood()
# )









