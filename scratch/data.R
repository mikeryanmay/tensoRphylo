# read tree
tree <- read.nexus("inst/testdata/sampled_ancestor_tree.nex")

# read data
data <- readNexusData("inst/testdata/sampled_ancestor_data.nex")

# add some nonsense
data <- data[-1,]
data <- rbind(data, "s_00" = c(0,1))

tp <- makeTensorPhylo(tree, data)
