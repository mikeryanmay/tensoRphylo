#' make the clado array function
CladogeneticEvents <- function(num_states) {
  return(new(CladoEvents, num_states))
}
