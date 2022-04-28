#' make the transition probability matrix function
TransitionMatrix <- function(num_states) {
  P <- new(ProbabilityMatrix, num_states)
  return(P)
}
