library(diversitree)

pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- NULL

while( is.null( phy ) ){
  phy <- tree.bisse(pars, max.t=30, x0=0, include.extinct=TRUE)
}

k.samples <- data.frame(taxon1="sp12", taxon2="sp12", timefrompresent=3.164384,
                        state=1, stringsAsFactors=FALSE)

hidden.states = FALSE
phy.k <- hisse:::AddKNodes(phy, k.samples)
fix.type <- hisse:::GetKSampleMRCA(phy.k, k.samples)
nb.tip <- Ntip(phy.k)
nb.node <- phy.k$Nnode
gen <- hisse:::FindGenerations(phy.k)

data <- data.frame(taxon=names(phy$tip.state), phy$tip.state,
                   stringsAsFactors=FALSE)
data <- hisse:::AddKData(data, k.samples)
data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
data.new <- data.new[phy.k$tip.label,]

dat.tab <- hisse:::OrganizeDataHiSSE(data.new, phy=phy.k, f=c(1,1),
                                     hidden.states=FALSE)
edge_details <- hisse:::GetEdgeDetails(phy=phy.k,
                                       intervening.intervals=strat.cache$intervening.intervals)
fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
pars.bisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01)

model.vec <- numeric(48)
model.vec[1:6] = pars.bisse
phy$node.label = NULL
cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states,
                                        nb.tip=Ntip(phy.k), nb.node=Nnode(phy.k), bad.likelihood=-300,
                                        f=c(1,1), ode.eps=0)

cache$psi <- 0.01
hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="herr_als",
                                    condition.on.survival=TRUE, root.p=c(1,1), node=fix.type$node,
                                    state=fix.type$state, fossil.taxa=fossil.taxa,
                                    fix.type=fix.type$type)
dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=1)
model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1,
                                      fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node,
                                      bad.likelihood=exp(-500), ode.eps=0)#
cache$psi <- 0.01
gen <- hisse:::FindGenerations(phy.k)
MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen,
                                    condition.on.survival=TRUE, root.type="herr_als", root.p=NULL,
                                    fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type)
# lambda: 0.10
# mu:     0.03
# eta:    0.01
# psi:    0.01

tp <- makeTensorPhylo(phy, nstates = 2)
tp$setApplyTreeLikCorrection(FALSE)
tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)

tp$setLambdaConstant(0.1)
tp$setMuConstant(0.03)
tp$setEtaConstantEqual(0.01)
tp$setPhiConstant(0.01)

tp$computeLogLikelihood() - MiSSE.logL + 2 * log(0.1)



