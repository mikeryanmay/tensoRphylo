#' make the clado array function
#' @internal
CladogeneticEvents <- function(num_states) {
  return(new(CladoEvents, num_states))
}
