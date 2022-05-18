#' Cladogenetic Event array
#'
#' The `CladoEvents` class.
#'
#' @details
#' This class is for sparse arrays of cladogenetic events.
#' The _i,j,k_'th element contains the probability that an ancestral lineage in state _i_ leaves one descendant in state _j_ and the other descendant in state _k_ during a speciation event.
#'
#' Objects of this class have only need to have a subset of possible events specified; unspecified events are assumed to occur with probability 0.
#' Because the elements of this array represent probabilities, the sum of all probabilities for a given ancestral state must be 1.
#' The class enforces that all values are between 0 and 1, and that ancestral probabilities sum to 1.
#'
#' For models with time-heterogeneous cladogenetic events, use the generic `c` to create a vector of event array (of class `CladoEventsList`):
#' \preformatted{
#' O_1 <- makeCladogeneticEvents(3)
#' O_2 <- makeCladogeneticEvents(3)
#' Os  <- c(O_1, O_2)
#' }
#'
#' ## Constructors
#' Create a new object of class `CladoEvents`.
#' - `makeCladogeneticEvents(dim)` where `dim` is an integer. Creates an empty array.
#' - `new(CladoEvents, dim)`. Equivalent to `makeCladogeneticEvents(dim)`.
#'
#' ## Accessors
#' `CladoEvents` provides intuitive access to elements of the array, but enforces some behavior.
#' - `[i,j,k]` returns value corresponding to an ancestor in state _i_ leaving one descendant in state_j_ and the other in states _k_.
#' - `[i,j,k] <- y` sets the value of the _ijk_th element to y. **_Will result in an error if you attempt to set a non-event (i = j = k), or attempt to set an element to a value outside of [0,1] bounds._**
#'
#' @examples
#' # create a cladogenetic event array.
#' W <- makeCladogeneticEvents(3)
#'
#' # set the probability of the ancestor in state 1 leaving left daughter in state 1 and right daughter in state 2
#' W[1,1,2] <- 0.2
#'
#' \dontrun{
#' # setting a diagonal value is prohibited
#' W[1,1,1] <- 0.1
#'
#' # likewise, setting a value outside of [0,1] is prohibited.
#' W[1,1,2] <- -0.1
#' W[1,1,2] <-  1.1
#' }
#'
#' # construct a time-homogeneous cladogenetic model
#' W  <- makeCladogeneticEvents(4)
#' tp <- new(TensorPhyloInstance, 4)
#' tp$setOmegaConstant(W)
#'
#' # construct a time-heterogeneous cladogenetic model
#' W_1 <- makeCladogeneticEvents(4)
#' W_2 <- makeCladogeneticEvents(4)
#' Ws  <- c(W_1, W_2)
#' t   <- 0.1
#' tp$setOmegaTimeDependent(t, Ws)
#' @name CladoEvents
#' @export
NULL

#' @export
makeCladogeneticEvents <- function(num_states) {
  return(new(CladoEvents, num_states))
}
