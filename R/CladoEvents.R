#' Cladogenetic Event array
#'
#' The `CladoEvents` class.
#'
#' @details
#' This class is for sparse arrays of cladogenetic events.
#' The _i,j,k_'th element containts the probability that an ancestral lineage in state _i_ leaves a left daughter in state _j_ and a right daughter in state _k_ during a speciation event.
#'
#' Objects of this class have only need to have a subset of possible events specified; unspecified events are assumed to occur with probability 0.
#' Because the elements of this array represent probabilities, the sum of all probabilities for a given ancestral state must be 1.
#' The class enforces that all values are between 0 and 1, and that ancestral probabilities sum to 1.
#'
#' For models with time-heterogeneous cladogenetic events, use the generic `c.CladoEvents` to create a vector of event array (of class `CladoEventsList`):
#' \preformatted{
#' O_1 <- makeCladogeneticEvents(3)
#' O_2 <- makeCladogeneticEvents(3)
#' Os  <- c(O_1, O_2)
#' }
#'
#' ## Constructors
#' Create a new object of class `CladoEvents`. For a more natural interface, use [`makeCladogeneticEvents`].
#' - `new(CladoEvents, dim)` where `dim` is an integer. Creates an empty array. Equivalent to `makeCladogeneticEvents(dim)`.
#'
#' ## Accessors
#' `CladoEvents` provides intuitive access to elements of the array, but enforces some behavior.
#' - `[i,j,k]` returns value corresponding to an ancestor in state _i_ leaving left and right descendants in states _j_ and _k_, respectively.
#' - `[i,j,k] <- y` sets the value of the _ijk_th element to y. **_Will result in an error if you attempt to set a non-event (i = j = k), or attempt to set an element to a value outside of [0,1] bounds._**
#'
#' @seealso [`makeCladogeneticEvents`]
#'
#' @examples
#' # create a cladogenetic event array.
#' O <- makeCladogeneticEvents(3)
#'
#' # set the probability of the ancestor in state 1 leaving left daughter in state 1 and right daughter in state 2
#' O[1,1,2] <- 0.2
#' \dontrun{
#' # setting a diagonal value is prohibited
#' O[1,1,1] <- 0.1
#'
#' # likewise, setting a value outside of [0,1] is prohibited.
#' O[1,1,2] <- -0.1
#' O[1,1,2] <-  1.1
#' }
#'
#' # construct a time-homogeneous cladogenetic model
#' O  <- makeCladogeneticEvents(4)
#' tp <- new(TensorPhylo, 4)
#' tp$setOmegaConstant(O)
#'
#' # construct a time-heterogeneous cladogenetic model
#' O_1 <- makeCladogeneticEvents(4)
#' O_2 <- makeCladogeneticEvents(4)
#' Os  <- c(O_1, O_2)
#' t   <- 0.1
#' tp$setOmegaTimeVarying(t, Os)
#' #' @name CladoEvents
#' @export
NULL

#' Cladogenetic array constructor.
#'
#' Creates a cladogenetic event  array
#'
#' @details
#' Creates an default object of class [`CladoEvents`].
#'
#' @param num_states The number of states.
#'
#' @return A default object of class [`CladoEvents`].
#'
#' @seealso [`CladoEvents`]
#'
#' @examples
#' # create a cladogenetic event array.
#' O <- makeCladogeneticEvents(3)
#'
#' # set the probability of the ancestor in state 1 leaving left daughter in state 1 and right daughter in state 2
#' O[1,1,2] <- 0.2
#'
#' @export
makeCladogeneticEvents <- function(num_states) {
  return(new(CladoEvents, num_states))
}
