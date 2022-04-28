# export some extra stuff
evalqOnLoad({

  setMethod("[", "Rcpp_CladoEvents", function(x, i, j, k) {
    x$getEventOneIndexed(c(i-1,j-1,k-1))
  })

  setMethod("[<-", "Rcpp_CladoEvents", function(x, i, j, k, value) {
    x$setEventOneIndexed(c(i-1,j-1,k-1), value)
    x
  })

  setMethod("c", "Rcpp_CladoEvents", function(x, ...) {

    # collect the arguments
    args <- list(...)

    # filter out invalid things
    args <- args[sapply(args, class) == "Rcpp_CladoEvents"]

    # create the new clado list
    l <- new(CladoEventsList, 1)
    l$addCladoEvents(x)

    # add any addition ones
    for(arg in args) {
      l$addCladoEvents(arg)
    }

    return(l)

  })

  setMethod("[", "Rcpp_RateMatrix", function(x, i, j) {
    x$getRateOneIndexed(i-1, j-1)
  })

  setMethod("[<-", "Rcpp_RateMatrix", function(x, i, j, value) {
    x$setRateOneIndexed(i-1, j-1, value)
    x
  })

  setMethod("c", "Rcpp_RateMatrix", function(x, ...) {

    # collect the arguments
    args <- list(...)

    # filter out invalid things
    args <- args[sapply(args, class) == "Rcpp_RateMatrix"]

    # create the new clado list
    l <- new(RateMatrixList, 1)
    l$addRateMatrix(x)

    # add any addition ones
    for(arg in args) {
      l$addRateMatrix(arg)
    }

    return(l)

  })

  setMethod("[", "Rcpp_ProbabilityMatrix", function(x, i, j) {
    x$getProbabilityOneIndexed(i-1, j-1)
  })

  setMethod("[<-", "Rcpp_ProbabilityMatrix", function(x, i, j, value) {
    x$setProbabilityOneIndexed(i-1, j-1, value)
    x
  })

  setMethod("c", "Rcpp_ProbabilityMatrix", function(x, ...) {

    # collect the arguments
    args <- list(...)

    # filter out invalid things
    args <- args[sapply(args, class) == "Rcpp_ProbabilityMatrix"]

    # create the new clado list
    l <- new(ProbabilityMatrixList, 1)
    l$addProbabilityMatrix(x)

    # add any addition ones
    for(arg in args) {
      l$addProbabilityMatrix(arg)
    }

    return(l)

  })

})
