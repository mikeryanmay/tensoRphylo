library(diversitree)

# simulate a tree with character data
pars <- c(0.1, 0.1, 0.03, 0.03, 0.05, 0.05)
repeat {

  # attempt to simulate a tree
  phy <- tree.bisse(pars, max.taxa = 20, x0 = 0, include.extinct = TRUE)

  # make sure tree is not null
  if ( is.null(phy) == FALSE ) {

    # make sure there is variation in the character
    if (length(unique(phy$tip.state)) == 2) {
      break
    }

  }

}

# add a root edge
phy$root.edge <- 0.5

# write the tree
write.nexus(phy, file="data-raw/tree.nex")

# write the character data
write.table(cbind(phy$tip.label, phy$tip.state), file = "data-raw/data.tsv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# read in the tree
test_tree <- read.nexus("data-raw/tree.nex")

# read in the data
test_data <- readDelimitedData("data-raw/data.tsv", delim = "\t", nstates = 2)

# drop the fossils
full_tree   <- test_tree
extant_tree <- drop.extinct(full_tree)
extant_tree$root.edge <- NULL

full_data   <- test_data
extant_data <- full_data[extant_tree$tip.label,]

# save the tree and data as RDA
save(full_tree, extant_tree, full_data, extant_data, file = "R/sysdata.Rda")





