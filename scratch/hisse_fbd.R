library(hisse)

# simulate...
phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.3, mu = 0.2)[[1]]
f <- GetFossils(phy, psi=0.05)
fbd.tree <- ProcessSimSample(phy, f)

# ugh. not sure I want to get this format into a proper newick string...
