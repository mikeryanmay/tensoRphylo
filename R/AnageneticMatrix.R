#' make the anagenetic rate matrix function
AnageneticMatrix <- function(num_states, rate = NULL) {
  if ( is.null(rate) == FALSE ) {
    Q <- new(RateMatrix, num_states, rate)
  } else {
    Q <- new(RateMatrix, num_states)
  }
  return(Q)
}
