library(diversitree)

# simulate tree
pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
names(pars) <- diversitree:::default.argnames.geosse()

set.seed(5)
phy <- tree.geosse(pars, max.t = 4, x0 = 0)

lik <- make.geosse(phy, phy$tip.state)
lik(pars, condition.surv = FALSE, root = ROOT.FLAT)

# statecols <- c("AB"="purple", "A"="blue", "B"="red")
# plot(phy, tip.color=statecols[phy$tip.state+1], cex=0.5)

# xA:  anagenetic change from A to AB
# xB:  anagenetic change from B to AB
# dA:  anagenetic change from AB to B
# dB:  anagenetic change from AB to A
# sA:  speciation in 1, no state change
# sB:  speciation in 2, no state change
# sAB: specation in 3, no state change

# 3->1,3 is sA
# 3->2,3 is sB
# 3->3,3 is sAB

# now set up in TP
map <- c("0" = "AB", "1" = "A", "2" = "B")

data <- matrix(0, length(phy$tip.state) , ncol = 3)
colnames(data) <- c("A","B","AB")
rownames(data) <- phy$tip.label
for(i in 1:length(phy$tip.state)) {
  data[i, map[as.character(phy$tip.state[i])]] <- 1.0
}

sA  <- pars[1]
sB  <- pars[2]
sAB <- pars[3]
xA  <- pars[4]
xB  <- pars[5]
dA  <- pars[6]
dB  <- pars[7]

# speciation rates per state
lambda <- c(sA, sB, sA + sB + sAB)

# extinction rates per state
mu <- c(dA, dB, 0)

# rate matrix
Q <- makeRateMatrix(3)
Q[3,1] <- xB
Q[3,2] <- xA
Q[1,3] <- dA
Q[2,3] <- dB

# cladogenetic events
O <- makeCladogeneticEvents(3)
O[3,3,1] <- 0.5 * sA  / (sA + sB + sAB)
O[3,1,3] <- 0.5 * sA  / (sA + sB + sAB)
O[3,3,2] <- 0.5 * sB  / (sA + sB + sAB)
O[3,2,3] <- 0.5 * sB  / (sA + sB + sAB)
O[3,1,2] <- 0.5 * sAB / (sA + sB + sAB)
O[3,2,1] <- 0.5 * sAB / (sA + sB + sAB)

tp <- makeTensorPhylo(phy, data)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$TIME)

tp$setLambdaStateVarying(lambda)
tp$setMuStateVarying(mu)
tp$setEtaConstantUnequal(Q)
tp$setOmegaConstant(O)

tp$computeLogLikelihood()

lik(pars, condition.surv = FALSE, root = ROOT.FLAT)


# do some manual correction
geo <- lik(pars, condition.surv = FALSE, root = ROOT.FLAT, intermediates = TRUE)
res <- attr(geo,'intermediates')
vals <- res$vals
lq <- res$lq
d.root <- vals[4:6]

sum(log(d.root)) + sum(lq)
as.numeric(geo)

diversitree:::rootfunc.geosse
diversitree:::rootfunc.classe








