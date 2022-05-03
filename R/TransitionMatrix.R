#' Probability Matrix
#'
#' The `ProbabilityMatrix` class.
#'
#' @details
#' This class is for transition-probability matrices, which describe the probability of changing from state _i_ (rows) to state _j_ (columns).
#' These matrices are useful for "synchronous" state-change events (where all lineages change simultaneously) at some time _t_.
#'
#' Objects of the class `ProbabilityMatrix` have elements that range from 0 to 1 such that the sum of each row is one.
#' Attempting to set the value of a diagonal element is prohibited, as they are computed automatically from the off-diagonal elements.
#' Likewise, the class enforces that element values `x` must be `0 <= x <= 1`.
#'
#' For models with time-heterogeneous rate matrices, use the generic `ProbabilityMatrix` to create a vector of rate matrices (of class `ProbabilityMatrixList`):
#' \preformatted{
#' P_1 <- makeProbabilityMatrix(3)
#' P_2 <- makeProbabilityMatrix(3)
#' Ps  <- c(P_1, P_2)
#' }
#'
#' ## Constructors
#' Create a new object of class `ProbabilityMatrix`.
#' - `makeProbabilityMatrix(dim)` where `dim` is an integer. Creates a `dim` x `dim` probability matrix with 1s along the diagonal.
#' - `new(ProbabilityMatrix, dim)`. Equivalent to `makeProbabilityMatrix(dim)`.
#'
#' ## Methods
#' Call a method on a `makeProbabilityMatrix` object `P` with `P$methodName(arguments)`.
#' - `getMatrix()`: Returns a standard R [base::matrix].
#'
#' ## Accessors
#' `ProbabilityMatrix` provides intuitive access to elements of the rate matrix, but enforces some behavior.
#' - `[i,j]` returns value in the _i_th row and _j_th column.
#' - `[i,j] <- y` sets the value of the _ij_th element to y. **_Will result in an error if you attempt to set a diagonal value, or attempt to set an off-diagonal element to below 0.0 or above 1.0._**
#'
#' @name ProbabilityMatrix
#'
#' @examples
#' # create a 4x4 probability matrix
#' P <- makeProbabilityMatrix(4)
#'
#' # set the rate from 1 to 2 to 0.2
#' P[1,2] <- 0.2
#'
#' \dontrun{
#' # setting a diagonal value is prohibited
#' P[1,1] <- 0.1
#'
#' # likewise, setting a value outside of [0,1] is prohibited.
#' P[1,2] <- -0.1
#' P[1,2] <- 1.1
#' }
#'
#' # construct a state-change model with one event at time t
#' # (we have to use a vector of P because there may be more than one
#' # time point, i.e., length(t) > 1).
#' P  <- makeProbabilityMatrix(4)
#' P[1,2] <- 0.1
#' t <- 0.1
#' tp <- new(TensorPhylo, 4)
#' tp$setZeta(t, c(P))
#'
#' @export
NULL

#' @export
makeProbabilityMatrix <- function(num_states) {
  P <- new(ProbabilityMatrix, num_states)
  return(P)
}
