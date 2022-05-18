library(tensoRphylo)
library(microbenchmark)

# read tree data
phy  <- read.nexus("scratch/tree.nex")
data <- readNexusData("scratch/data.nex")
data[,2] <- 1 - data[,1]
data <- data[,1:2]

# parameters
lambda <- 1
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

# tp$computeLogLikelihood()
# tp$setOmegaConstant(O)
# tp$computeLogLikelihood()

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

# make the classe model
data_vec  <- (data %*% c(0,1))[,1] + 1
model <- diversitree::make.classe(phy, data_vec, k=2, sampling.f=c(1,1), control = list(tol = 1e-8))

sprintf("%.10f", tp$computeLogLikelihood())
sprintf("%.10f", model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT))

microbenchmark(
  tp$computeLogLikelihood(),
  model(pars, condition.surv = FALSE, root = diversitree::ROOT.FLAT),
  times = 100L
)














# path to tree file
# tree_file <- system.file("testdata", "extant_tree.nex",
#                          package = "tensoRphylo")
#
# # read the tree
# tree <- ape::read.nexus(tree_file)
#
# # path to data file
# data_file <- system.file("testdata", "extant_data.nex",
#                          package = "tensoRphylo")
#
# # read the data
# data <- readNexusData(data_file)

tp <- makeTensorPhylo(phy, data)
tp$setLambdaConstant(0.5)
tp$setMuConstant(0.000001)
tp$setEtaConstantEqual(0.1)

# make an empty cladogenetic array.
O <- makeCladogeneticEvents(num_states = 2)

# print the empty array
O

# specify a few cladogenetic events
O[1,1,2] <- 0.2
O[2,2,1] <- 0.3

# print the populated array
O

# compute the likelihood
tp$computeLogLikelihood()

# set the cladogenetic array
tp$setOmegaConstant(O)

# compute the likelihood
tp$computeLogLikelihood()
