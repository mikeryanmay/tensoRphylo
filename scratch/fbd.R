source("scratch/fbd_equations.R")
source("scratch/bdss_prob.R")

# read tree
tree <- tensoRphylo:::full_tree

# parameters
lambda <- 0.1
mu     <- 0.03
phi    <- 0.05
rho    <- 1.0

# compute prob
fbd_likelihood_stem(tree, lambda, mu, phi, rho)

# compute with TP
tp <- makeTensorPhylo(tree, nstates = 2)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType( conditionalProbability$STEM_ONE_SAMPLE )

# set parameters
tp$setLambdaConstant(lambda)
tp$setMuConstant(mu)
tp$setDeltaConstant(phi)

# compute likelihood
tp$computeLogLikelihood()

# try J's function
ages         <- tree.age(tree, digits = 6)
branch_times <- ages[-c(1:length(tree$tip.label)),1]
branch_times <- c(branch_times, max(branch_times) + phy$root.edge)
sample_times <- ages$ages[ages$ages > 0]
sample_times <- sample_times[sample_times %in% branch_times == FALSE]

tree_logprob(branch_times, sample_times, root = 0, survival = TRUE, lambda, mu, phi, rho = 1)

