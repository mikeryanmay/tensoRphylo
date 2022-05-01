library(tensoRphylo)
library(microbenchmark)
library(RevGadgets)
library(phytools)

# read some stuff
tree <- read.nexus("test/tree.nex")
data <- readNexusData("test/data.nex")

# create the tp object
tp   <- tensoRphylo(tree, data)

tp$setLambdaConstant(0.01)
tp$setEtaConstantEqual(0.01)
tp$computeLogLikelihood()

# simulate a history
histories <- tp$drawStochasticMaps(10)

# plot ancestral-state estimates
ss <- summary(histories)
plot(ss)
