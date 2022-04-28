#' make the anagenetic rate matrix function
#' @internal
AnageneticMatrix <- function(num_states, rate = NULL) {
  if ( is.null(rate)  ) {
    Q <- new(RateMatrix, num_states)
  } else {
    Q <- new(RateMatrix, num_states, rate)
  }
  return(Q)
}
