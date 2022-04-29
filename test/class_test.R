library(microbenchmark)

# read some stuff
tree <- read.nexus("test/tree.nex")
data <- readNexusData("test/data.nex")
tp   <- tensoRphylo(tree, data)

tp$setLambdaConstant(0.01)
tp$setEtaConstantEqual(0.01)
tp$computeLogLikelihood()

history <- tp$drawStochasticMap()

