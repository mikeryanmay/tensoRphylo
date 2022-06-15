library(ape)
library(FossilSim)

# simulate a tree with zero-length branches
repeat {
  tree <- sim.fbd.age(4, 1, 1, 0.5, 0.5)[[1]]
  if ( inherits(tree, "phylo") ) {
    break
  }
}

tree <- ladderize(tree)

# convert it to a SA tree (degree-one nodes)
new_tree <- tensoRphylo:::.zero.bl.to.SA.tree(tree)

# plot the sampled ancestors with names
par(mfrow = c(2,1))
plot(tree, no.margin = TRUE, root.edge = TRUE)

plot(new_tree, no.margin = TRUE, root.edge = TRUE)
sampled_ancestors <- (1:new_tree$Nnode + Ntip(new_tree))[new_tree$node.label != ""]
nodelabels(node = sampled_ancestors, text = new_tree$node.label[new_tree$node.label != ""], adj = c(-0.1,0), frame = "none", srt = 45)
nodelabels(node = sampled_ancestors, pch  = 19)













