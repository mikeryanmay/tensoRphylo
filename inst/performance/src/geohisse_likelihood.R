library(hisse)

GeoHiSSELikelihood <- function(phy, data, f=c(1,1,1), turnover=c(1,2,3), eps=c(1,2), hidden.states=FALSE, trans.rate=NULL, assume.cladogenetic=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=TRUE, sann.its=1000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, mag.san.start=0.5, starting.vals=NULL, turnover.upper=1000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, dt.threads=1){

  ## Temporary fix for the current BUG:
  if( !is.null(phy$node.label) ) phy$node.label <- NULL

  if(!is.null(root.p)) {
    ## The vector of the root.p need to be as long as the speciation vector.
    if(hidden.states == TRUE){
      if( length( root.p ) == 3 ){
        root.p <- rep(root.p, 10)
        root.p <- root.p / sum(root.p)
        warning("For hidden areas, you need to specify the root.p for all hidden areas. We have adjusted it so that there's equal chance for among all hidden areas.")
      } else{
        root.p.new <- numeric(30)
        root.p.new[1:length(root.p)] <- root.p
        root.p <- root.p.new
        root.p <- root.p / sum(root.p)
      }
    }else{
      ## All good:
      root.p <- root.p / sum(root.p)
    }

  }

  data.table::setDTthreads(threads=dt.threads)

  if(sann == FALSE & is.null(starting.vals)){
    warning("You have chosen to rely on the internal starting points that generally work but does not guarantee finding the MLE.")
  }

  if(!root.type == "madfitz" & !root.type == "herr_als"){
    stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
  }

  if(is.null(trans.rate)){
    stop("Rate matrix needed. See TransMatMakerGeoHiSSE() to create one.")
  }

  if(hidden.states == TRUE & dim(trans.rate)[1]<3){
    stop("You chose a hidden state but this is not reflected in the transition matrix")
  }

  ## Return error message if the data is not in the correct format.
  if( !inherits(data, what = c("matrix","data.frame")) ){
    stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
  }
  if( !ncol( data ) == 2 ){
    stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
  }
  ## Check if the states are %in% 0:2:
  states.check <- all( as.numeric(data[,2]) %in% 0:2 )
  if( !states.check ){
    stop("states need to be one of 0, 1, or 2. See help.")
  }

  ## Check if 'hidden.states' parameter is congruent with the turnover vector:
  if( length(turnover) > 3 & !hidden.states ){
    stop("Turnover has more than 3 elements but 'hidden.states' was set to FALSE. Please set 'hidden.states' to TRUE if the model include more than one rate class.")
  }
  if( length(turnover) == 3 & hidden.states ){
    stop("Turnover has only 3 elements but 'hidden.states' was set to TRUE. Please set 'hidden.states' to FALSE if the model does not include hidden rate classes.")
  }

  pars <- numeric(380)

  if(dim(trans.rate)[2]==3){
    rate.cats <- 1
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    trans.tmp <- c(trans.rate["(00)", "(11)"], trans.rate["(00)", "(01)"], trans.rate["(11)", "(00)"], trans.rate["(11)", "(01)"],  trans.rate["(01)", "(00)"],  trans.rate["(01)", "(11)"])
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    category.rates.unique <- 0
    pars.tmp <- c(pars.tmp, trans.tmp)
    pars[1:11] <- pars.tmp
  }

  if(dim(trans.rate)[2]==6){
    rate.cats <- 2
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 1), rep("(11A)", 1), rep("(01A)", 1), rep("(00B)", 1), rep("(11B)", 1), rep("(01B)", 1))
    cols <- c("(00B)", "(11B)", "(01B)", "(00A)", "(11A)", "(01A)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1], rep(0,8), category.rate.shift[2], rep(0,8), category.rate.shift[3], rep(0,8))
    category.rate.shiftB <- c(category.rate.shift[4], rep(0,8), category.rate.shift[5], rep(0,8), category.rate.shift[6], rep(0,8))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB)
    pars[1:length(pars.tmp)] <- pars.tmp
  }

  if(dim(trans.rate)[2]==9){
    rate.cats <- 3
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 2), rep("(11A)", 2), rep("(01A)", 2), rep("(00B)", 2), rep("(11B)", 2), rep("(01B)", 2), rep("(00C)", 2), rep("(11C)", 2), rep("(01C)", 2))
    cols <- c("(00B)", "(00C)", "(11B)", "(11C)", "(01B)", "(01C)", "(00A)", "(00C)", "(11A)", "(11C)", "(01A)", "(01C)", "(00A)", "(00B)", "(11A)", "(11B)", "(01A)", "(01B)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:2], rep(0,7), category.rate.shift[3:4], rep(0,7), category.rate.shift[5:6], rep(0,7))
    category.rate.shiftB <- c(category.rate.shift[7:8], rep(0,7), category.rate.shift[9:10], rep(0,7), category.rate.shift[11:12], rep(0,7))
    category.rate.shiftC <- c(category.rate.shift[13:14], rep(0,7), category.rate.shift[15:16], rep(0,7), category.rate.shift[17:18], rep(0,7))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC)
    pars[1:length(pars.tmp)] <- pars.tmp
  }

  if(dim(trans.rate)[2]==12){
    rate.cats <- 4
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 3), rep("(11A)", 3), rep("(01A)", 3), rep("(00B)", 3), rep("(11B)", 3), rep("(01B)", 3), rep("(00C)", 3), rep("(11C)", 3), rep("(01C)", 3), rep("(00D)", 3), rep("(11D)", 3), rep("(01D)", 3))
    cols <- c("(00B)", "(00C)", "(00D)", "(11B)", "(11C)", "(11D)", "(01B)", "(01C)", "(01D)", "(00A)", "(00C)", "(00D)", "(11A)", "(11C)", "(11D)", "(01A)", "(01C)", "(01D)", "(00A)", "(00B)", "(00D)", "(11A)", "(11B)", "(11D)", "(01A)", "(01B)", "(01D)", "(00A)", "(00B)", "(00C)", "(11A)", "(11B)", "(11C)", "(01A)", "(01B)", "(01C)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:3], rep(0,6), category.rate.shift[4:6], rep(0,6), category.rate.shift[7:9], rep(0,6))
    category.rate.shiftB <- c(category.rate.shift[10:12], rep(0,6), category.rate.shift[13:15], rep(0,6), category.rate.shift[16:18], rep(0,6))
    category.rate.shiftC <- c(category.rate.shift[19:21], rep(0,6), category.rate.shift[22:24], rep(0,6), category.rate.shift[25:27], rep(0,6))
    category.rate.shiftD <- c(category.rate.shift[28:30], rep(0,6), category.rate.shift[31:33], rep(0,6), category.rate.shift[34:36], rep(0,6))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD)
    pars[1:length(pars.tmp)] <- pars.tmp
  }

  if(dim(trans.rate)[2]==15){
    rate.cats <- 5
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 4), rep("(11A)", 4), rep("(01A)", 4), rep("(00B)", 4), rep("(11B)", 4), rep("(01B)", 4), rep("(00C)", 4), rep("(11C)", 4), rep("(01C)", 4), rep("(00D)", 4), rep("(11D)", 4), rep("(01D)", 4), rep("(00E)", 4), rep("(11E)", 4), rep("(01E)", 4))
    cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(11B)", "(11C)", "(11D)", "(11E)", "(01B)", "(01C)", "(01D)", "(01E)", "(00A)", "(00C)", "(00D)", "(00E)", "(11A)", "(11C)", "(11D)", "(11E)", "(01A)", "(01C)", "(01D)", "(01E)", "(00A)", "(00B)", "(00D)", "(00E)", "(11A)", "(11B)", "(11D)", "(11E)", "(01A)", "(01B)", "(01D)", "(01E)", "(00A)", "(00B)", "(00C)", "(00E)", "(11A)", "(11B)", "(11C)", "(11E)", "(01A)", "(01B)", "(01C)", "(01E)", "(00A)", "(00B)", "(00C)", "(00D)", "(11A)", "(11B)", "(11C)", "(11D)", "(01A)", "(01B)", "(01C)", "(01D)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:4], rep(0,5), category.rate.shift[5:8], rep(0,5), category.rate.shift[9:12], rep(0,5))
    category.rate.shiftB <- c(category.rate.shift[13:16], rep(0,5), category.rate.shift[17:20], rep(0,5), category.rate.shift[21:24], rep(0,5))
    category.rate.shiftC <- c(category.rate.shift[25:28], rep(0,5), category.rate.shift[29:32], rep(0,5), category.rate.shift[33:36], rep(0,5))
    category.rate.shiftD <- c(category.rate.shift[37:40], rep(0,5), category.rate.shift[41:44], rep(0,5), category.rate.shift[45:48], rep(0,5))
    category.rate.shiftE <- c(category.rate.shift[49:52], rep(0,5), category.rate.shift[53:56], rep(0,5), category.rate.shift[57:60], rep(0,5))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE)
    pars[1:length(pars.tmp)] <- pars.tmp
  }

  if(dim(trans.rate)[2]==18){
    rate.cats <- 6
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 5), rep("(11A)", 5), rep("(01A)", 5), rep("(00B)", 5), rep("(11B)", 5), rep("(01B)", 5), rep("(00C)", 5), rep("(11C)", 5), rep("(01C)", 5), rep("(00D)", 5), rep("(11D)", 5), rep("(01D)", 5), rep("(00E)", 5), rep("(11E)", 5), rep("(01E)", 5), rep("(00F)", 5), rep("(11F)", 5), rep("(01F)", 5))
    cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:5], rep(0,4), category.rate.shift[6:10], rep(0,4), category.rate.shift[11:15], rep(0,4))
    category.rate.shiftB <- c(category.rate.shift[16:20], rep(0,4), category.rate.shift[21:25], rep(0,4), category.rate.shift[26:30], rep(0,4))
    category.rate.shiftC <- c(category.rate.shift[31:35], rep(0,4), category.rate.shift[36:40], rep(0,4), category.rate.shift[41:45], rep(0,4))
    category.rate.shiftD <- c(category.rate.shift[46:50], rep(0,4), category.rate.shift[51:55], rep(0,4), category.rate.shift[56:60], rep(0,4))
    category.rate.shiftE <- c(category.rate.shift[61:65], rep(0,4), category.rate.shift[66:70], rep(0,4), category.rate.shift[71:75], rep(0,4))
    category.rate.shiftF <- c(category.rate.shift[76:80], rep(0,4), category.rate.shift[81:85], rep(0,4), category.rate.shift[86:90], rep(0,4))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF)
    pars[1:length(pars.tmp)] <- pars.tmp
  }


  if(dim(trans.rate)[2]==21){
    rate.cats <- 7
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 6), rep("(11A)", 6), rep("(01A)", 6), rep("(00B)", 6), rep("(11B)", 6), rep("(01B)", 6), rep("(00C)", 6), rep("(11C)", 6), rep("(01C)", 6), rep("(00D)", 6), rep("(11D)", 6), rep("(01D)", 6), rep("(00E)", 6), rep("(11E)", 6), rep("(01E)", 6), rep("(00F)", 6), rep("(11F)", 6), rep("(01F)", 6), rep("(00G)", 6), rep("(11G)", 6), rep("(01G)", 6))
    cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:6], rep(0,3), category.rate.shift[7:12], rep(0,3), category.rate.shift[13:18], rep(0,3))
    category.rate.shiftB <- c(category.rate.shift[19:24], rep(0,3), category.rate.shift[25:30], rep(0,3), category.rate.shift[31:36], rep(0,3))
    category.rate.shiftC <- c(category.rate.shift[37:42], rep(0,3), category.rate.shift[43:48], rep(0,3), category.rate.shift[49:54], rep(0,3))
    category.rate.shiftD <- c(category.rate.shift[55:60], rep(0,3), category.rate.shift[61:66], rep(0,3), category.rate.shift[67:72], rep(0,3))
    category.rate.shiftE <- c(category.rate.shift[73:78], rep(0,3), category.rate.shift[79:84], rep(0,3), category.rate.shift[85:90], rep(0,3))
    category.rate.shiftF <- c(category.rate.shift[91:96], rep(0,3), category.rate.shift[97:102], rep(0,3), category.rate.shift[103:108], rep(0,3))
    category.rate.shiftG <- c(category.rate.shift[109:114], rep(0,3), category.rate.shift[115:120], rep(0,3), category.rate.shift[121:126], rep(0,3))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG)
    pars[1:length(pars.tmp)] <- pars.tmp
  }


  if(dim(trans.rate)[2]==24){
    rate.cats <- 8
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)", "(00H)", "(00H)",  "(11H)", "(11H)",  "(01H)", "(01H)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)",  "(11H)", "(01H)", "(00H)", "(01H)", "(00H)",  "(11H)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 7), rep("(11A)", 7), rep("(01A)", 7), rep("(00B)", 7), rep("(11B)", 7), rep("(01B)", 7), rep("(00C)", 7), rep("(11C)", 7), rep("(01C)", 7), rep("(00D)", 7), rep("(11D)", 7), rep("(01D)", 7), rep("(00E)", 7), rep("(11E)", 7), rep("(01E)", 7), rep("(00F)", 7), rep("(11F)", 7), rep("(01F)", 7), rep("(00G)", 7), rep("(11G)", 7), rep("(01G)", 7), rep("(00H)", 7), rep("(11H)", 7), rep("(01H)", 7))
    cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:7], rep(0,2), category.rate.shift[8:14], rep(0,2), category.rate.shift[15:21], rep(0,2))
    category.rate.shiftB <- c(category.rate.shift[22:28], rep(0,2), category.rate.shift[29:35], rep(0,2), category.rate.shift[36:42], rep(0,2))
    category.rate.shiftC <- c(category.rate.shift[43:49], rep(0,2), category.rate.shift[50:56], rep(0,2), category.rate.shift[57:63], rep(0,2))
    category.rate.shiftD <- c(category.rate.shift[64:70], rep(0,2), category.rate.shift[71:77], rep(0,2), category.rate.shift[78:84], rep(0,2))
    category.rate.shiftE <- c(category.rate.shift[85:91], rep(0,2), category.rate.shift[92:98], rep(0,2), category.rate.shift[99:105], rep(0,2))
    category.rate.shiftF <- c(category.rate.shift[106:112], rep(0,2), category.rate.shift[113:119], rep(0,2), category.rate.shift[120:126], rep(0,2))
    category.rate.shiftG <- c(category.rate.shift[127:133], rep(0,2), category.rate.shift[134:140], rep(0,2), category.rate.shift[141:147], rep(0,2))
    category.rate.shiftH <- c(category.rate.shift[148:154], rep(0,2), category.rate.shift[155:161], rep(0,2), category.rate.shift[162:168], rep(0,2))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG, category.rate.shiftH)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG, turnover[22:24], extinct.frac.tmp[15:16], trans.tmp[43:48], category.rate.shiftH)
    pars[1:length(pars.tmp)] <- pars.tmp
  }


  if(dim(trans.rate)[2]==27){
    rate.cats <- 9
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)", "(00H)", "(00H)",  "(11H)", "(11H)",  "(01H)", "(01H)", "(00I)", "(00I)",  "(11I)", "(11I)",  "(01I)", "(01I)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)",  "(11H)", "(01H)", "(00H)", "(01H)", "(00H)",  "(11H)",  "(11I)", "(01I)", "(00I)", "(01I)", "(00I)",  "(11I)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 8), rep("(11A)", 8), rep("(01A)", 8), rep("(00B)", 8), rep("(11B)", 8), rep("(01B)", 8), rep("(00C)", 8), rep("(11C)", 8), rep("(01C)", 8), rep("(00D)", 8), rep("(11D)", 8), rep("(01D)", 8), rep("(00E)", 8), rep("(11E)", 8), rep("(01E)", 8), rep("(00F)", 8), rep("(11F)", 8), rep("(01F)", 8), rep("(00G)", 8), rep("(11G)", 8), rep("(01G)", 8), rep("(00H)", 8), rep("(11H)", 8), rep("(01H)", 8), rep("(00I)", 8), rep("(11I)", 8), rep("(01I)", 8))
    cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- c(category.rate.shift[1:8], rep(0,1), category.rate.shift[9:16], rep(0,1), category.rate.shift[17:24], rep(0,1))
    category.rate.shiftB <- c(category.rate.shift[25:32], rep(0,1), category.rate.shift[33:40], rep(0,1), category.rate.shift[41:48], rep(0,1))
    category.rate.shiftC <- c(category.rate.shift[49:56], rep(0,1), category.rate.shift[57:64], rep(0,1), category.rate.shift[65:72], rep(0,1))
    category.rate.shiftD <- c(category.rate.shift[73:80], rep(0,1), category.rate.shift[81:88], rep(0,1), category.rate.shift[89:96], rep(0,1))
    category.rate.shiftE <- c(category.rate.shift[97:104], rep(0,1), category.rate.shift[105:112], rep(0,1), category.rate.shift[113:120], rep(0,1))
    category.rate.shiftF <- c(category.rate.shift[121:128], rep(0,1), category.rate.shift[129:136], rep(0,1), category.rate.shift[137:144], rep(0,1))
    category.rate.shiftG <- c(category.rate.shift[145:152], rep(0,1), category.rate.shift[153:160], rep(0,1), category.rate.shift[161:168], rep(0,1))
    category.rate.shiftH <- c(category.rate.shift[169:176], rep(0,1), category.rate.shift[177:184], rep(0,1), category.rate.shift[185:192], rep(0,1))
    category.rate.shiftI <- c(category.rate.shift[193:200], rep(0,1), category.rate.shift[201:208], rep(0,1), category.rate.shift[209:216], rep(0,1))
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG, category.rate.shiftH, category.rate.shiftI)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG, turnover[22:24], extinct.frac.tmp[15:16], trans.tmp[43:48], category.rate.shiftH, turnover[25:27], extinct.frac.tmp[17:18], trans.tmp[49:54], category.rate.shiftI)
    pars[1:length(pars.tmp)] <- pars.tmp
  }

  if(dim(trans.rate)[2]==30){
    rate.cats <- 10
    pars.tmp <- turnover
    extinct.frac.tmp <- eps
    extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, extinct.frac.tmp)
    for.late.adjust <- max(pars.tmp)
    rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)", "(00H)", "(00H)",  "(11H)", "(11H)",  "(01H)", "(01H)", "(00I)", "(00I)",  "(11I)", "(11I)",  "(01I)", "(01I)", "(00J)", "(00J)",  "(11J)", "(11J)",  "(01J)", "(01J)")
    cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)",  "(11H)", "(01H)", "(00H)", "(01H)", "(00H)",  "(11H)",  "(11I)", "(01I)", "(00I)", "(01I)", "(00I)",  "(11I)",  "(11J)", "(01J)", "(00J)", "(01J)", "(00J)", "(11J)")
    trans.tmp <- trans.rate[cbind(rows,cols)]
    trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, trans.tmp)
    rows <- c(rep("(00A)", 9), rep("(11A)", 9), rep("(01A)", 9), rep("(00B)", 9), rep("(11B)", 9), rep("(01B)", 9), rep("(00C)", 9), rep("(11C)", 9), rep("(01C)", 9), rep("(00D)", 9), rep("(11D)", 9), rep("(01D)", 9), rep("(00E)", 9), rep("(11E)", 9), rep("(01E)", 9), rep("(00F)", 9), rep("(11F)", 9), rep("(01F)", 9), rep("(00G)", 9), rep("(11G)", 9), rep("(01G)", 9), rep("(00H)", 9), rep("(11H)", 9), rep("(01H)", 9), rep("(00I)", 9), rep("(11I)", 9), rep("(01I)", 9), rep("(00J)", 9), rep("(11J)", 9), rep("(01J)", 9))
    cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01I)", "(01J)",  "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)")
    category.tmp <- trans.rate[cbind(rows,cols)]
    category.tmp[category.tmp==0] <- NA
    category.rate.shift <- category.tmp + for.late.adjust
    category.rate.shift[is.na(category.rate.shift)] <- 0
    category.rate.shiftA <- category.rate.shift[1:27]
    category.rate.shiftB <- category.rate.shift[28:54]
    category.rate.shiftC <- category.rate.shift[55:81]
    category.rate.shiftD <- category.rate.shift[82:108]
    category.rate.shiftE <- category.rate.shift[109:135]
    category.rate.shiftF <- category.rate.shift[136:162]
    category.rate.shiftG <- category.rate.shift[163:189]
    category.rate.shiftH <- category.rate.shift[190:216]
    category.rate.shiftI <- category.rate.shift[217:243]
    category.rate.shiftJ <- category.rate.shift[244:270]
    category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG, category.rate.shiftH, category.rate.shiftI, category.rate.shiftJ)
    category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
    pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG, turnover[22:24], extinct.frac.tmp[15:16], trans.tmp[43:48], category.rate.shiftH, turnover[25:27], extinct.frac.tmp[17:18], trans.tmp[49:54], category.rate.shiftI, turnover[28:30], extinct.frac.tmp[19:20], trans.tmp[55:60], category.rate.shiftJ)
    pars[1:length(pars.tmp)] <- pars.tmp
  }

  np <- max(pars)
  pars[pars==0] <- np+1

  data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
  data.new <- data.new[phy$tip.label,]

  #This is used to scale starting values to account for sampling:
  if(length(f) == 3){
    freqs <- table(data.new[,1])
    if(length(freqs == 2)){
      samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f[as.numeric(names(freqs))+1])
    }else{
      samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f[as.numeric(names(freqs))+1])
    }
  }else{
    if(length(f) == Ntip(phy)){
      stop("This is functionality has been temporarily removed.")
      #samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / mean(f[as.numeric(names(freqs))+1]))
    }else{
      stop("The vector of sampling frequencies does not match the number of tips in the tree.")
    }
  }

  if(is.null(restart.obj)){
    if(sum(eps)==0){
      init.pars <- hisse:::starting.point.geosse(phy, eps=0, samp.freq.tree=samp.freq.tree)
    }else{
      init.pars <- hisse:::starting.point.geosse(phy, eps=mag.san.start, samp.freq.tree=samp.freq.tree)
      #init.pars <- starting.point.generator(phy, k=3, samp.freq.tree=samp.freq.tree, q.div=5, yule=FALSE)
    }
    names(init.pars) <- NULL

    if(is.null(starting.vals)){
      #def.set.pars <- rep(c(log(init.pars[1:3]), log(init.pars[4:5]), rep(log(init.pars[6:7]*.1),3), rep(log(.01), 27)), rate.cats)
      #def.set.pars <- rep(c(log(init.pars[1:3]), log(init.pars[4:5]), rep(log(init.pars[7:8]*.1),3), rep(log(.01), 27)), rate.cats)
      def.set.pars <- rep(c(log(init.pars[1]+init.pars[4]), log(init.pars[2]+init.pars[5]), log(sum(init.pars[1:3])), log(init.pars[4]/init.pars[1]),  log(init.pars[5]/init.pars[2]), rep(log(init.pars[6:7]),3), rep(log(0.01), 27)), rate.cats)
    }else{
      ## Check the format for the starting.vals and accepts legacy mode if necessary.
      if( !length(starting.vals) %in% c(3,7) ){
        stop("Wrong length of starting.vals")
      }
      if( length(starting.vals) == 7 ){
        cat("Using developer mode for starting.vals.", "\n")
        def.set.pars <- rep(c(log(starting.vals[1:3]), log(starting.vals[4:5]), rep(log(starting.vals[6:7]),3), rep(log(0.01), 27)), rate.cats)
      } else{
        def.set.pars <- rep(c(log( rep(starting.vals[1],3) ), log( rep(starting.vals[2],2) ), rep(log(starting.vals[3]),6), rep(log(0.01), 27)), rate.cats)
      }
    }
    if(bounded.search == TRUE){
      upper.full <- rep(c(rep(log(turnover.upper),3), rep(log(eps.upper),2), rep(log(trans.upper),2*3), rep(log(10), 27)), rate.cats)
      #upper.full <- rep(c(rep(log(10000),3), rep(log(3),2), rep(log(trans.upper), 2*3), rep(log(10), 27)), rate.cats)
    }else{
      upper.full <- rep(21,length(def.set.pars))
    }
    np.sequence <- 1:np
    ip <- numeric(np)
    upper <- numeric(np)
    for(i in np.sequence){
      ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
      upper[i] <- upper.full[which(pars == np.sequence[i])[1]]
    }
    lower <- rep(-20, length(ip))
  }else{
    upper <- restart.obj$upper.bounds
    lower <- restart.obj$lower.bounds
    pars <- restart.obj$index.par
    ip <- numeric(length(unique(restart.obj$index.par))-1)
    for(k in 1:length(ip)){
      ip[k] <- log(restart.obj$solution[which(restart.obj$index.par==k)][1])
    }
  }

  # Some new prerequisites #
  gen <- hisse:::FindGenerations(phy)
  dat.tab <- hisse:::OrganizeDataGeo(data=data.new[,1], phy=phy, f=f, hidden.states=hidden.states)
  nb.tip <- Ntip(phy)
  nb.node <- phy$Nnode
  ##########################

  res <- function() {
    -hisse:::DevOptimizeGeoHiSSEfast(ip, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
  }

  return(res)

}


assignInNamespace("GeoHiSSE", GeoHiSSELikelihood, "hisse")
