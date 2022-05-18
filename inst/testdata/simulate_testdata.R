library(ape)
library(geiger)
library(tensoRphylo)

# read test tree
sa_tree <- read.nexus("sampled_ancestor_tree.nex")
sa_data <- readNexusData("sampled_ancestor_data.nex")

# prune tree to serial-sampled tree
ss_tree <- collapse.singles(sa_tree)
ss_data <- sa_data[ss_tree$tip.label,]

# prune tree to have only extant samples
extant_tree <- drop.extinct(ss_tree)
extant_tree$root.edge <- NULL
extant_data <- sa_data[extant_tree$tip.label,]

# create the data vector
data_vec <- (sa_data %*% c(0,1))[,1]
data_tmp <- as.character(data_vec)
names(data_tmp) <- names(data_vec)
data_tmp[1] <- "0/1"
data_vec <- data_tmp

extant_data_vec <- (extant_data %*% c(0,1))[,1]
extant_data_tmp <- as.character(extant_data_vec)
names(extant_data_tmp) <- names(extant_data_vec)
extant_data_tmp[1] <- "0/1"
extant_data_vec <- extant_data_tmp

# save the data
saveRDS(sa_tree, file         = "sampled_ancestor_tree.Rda")
saveRDS(sa_data, file         = "sampled_ancestor_data.Rda")
saveRDS(data_vec, file        = "sampled_ancestor_data_vec.Rda")
saveRDS(ss_tree, file         = "serial_sampled_tree.Rda")
saveRDS(ss_data, file         = "serial_sampled_data.Rda")
saveRDS(extant_tree, file     = "extant_tree.Rda")
saveRDS(extant_data, file     = "extant_data.Rda")
saveRDS(extant_data_vec, file = "extant_data_vec.Rda")
