library(ape)
library(geiger)
library(tensoRphylo)

# read test tree
sa_tree <- read.nexus("tests/testthat/testdata/sampled_ancestor_tree.nex")
sa_data <- readNexusData("tests/testthat/testdata/sampled_ancestor_data.nex")

# prune tree to serial-sampled tree
ss_tree <- collapse.singles(sa_tree)
ss_data <- sa_data[ss_tree$tip.label,]

# prune tree to have only extant samples
extant_tree <- drop.extinct(ss_tree)
extant_tree$root.edge <- NULL
extant_data <- sa_data[extant_tree$tip.label,]

# save the data
saveRDS(sa_tree, file = "tests/testthat/testdata/sampled_ancestor_tree.Rda")
saveRDS(sa_data, file = "tests/testthat/testdata/sampled_ancestor_data.Rda")
saveRDS(ss_tree, file = "tests/testthat/testdata/serial_sampled_tree.Rda")
saveRDS(ss_data, file = "tests/testthat/testdata/serial_sampled_data.Rda")
saveRDS(extant_tree, file = "tests/testthat/testdata/extant_tree.Rda")
saveRDS(extant_data, file = "tests/testthat/testdata/extant_data.Rda")

# readDelimitedData("tests/testthat/testdata/sampled_ancestor_data.tsv", delim = "\t", nstates = 2)
# readDelimitedData("tests/testthat/testdata/sampled_ancestor_data.csv", delim = ",", nstates = 2)
