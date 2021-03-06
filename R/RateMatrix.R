#' Rate Matrix
#'
#' The `RateMatrix` class.
#'
#' @details
#' This class is for instantaneous rate matrices, which describe the rate of change from state _i_ (rows) to state _j_ (columns).
#' In tensorphylo, these matrices are used to describe "anagenetic", "asynchronous" events, i.e., changes in state that are not associated with speciation events, and occur independently among lineages.
#' This is directly analogous to the Q matrix of a substitution model.
#'
#' Objects of the class `RateMatrix` have off-diagonal elements that are strictly positive and diagonal elements such that the sum of each row is zero.
#' Attempting to set the value of a diagonal element is prohibited, as they are computed automatically from the off-diagonal elements.
#' Likewise, the class enforces that off-diagonal elements must >= 0.
#'
#' For models with time-dependent rate matrices, use the generic `c` to create a vector of rate matrices (of class `RateMatrixList`):
#' \preformatted{
#' Q_1 <- makeRateMatrix(3, 0.1)
#' Q_2 <- makeRateMatrix(3, 0.2)
#' Qs  <- c(Q_1, Q_2)
#' }
#'
#' ## Constructors
#' Create a new object of class `RateMatrix`.
#' - `makeRateMatrix(num_states)` where `num_states` is an integer. Creates a `num_states` x `num_states` rate matrix with rate zero.
#' - `makeRateMatrix(num_states, rate)` where `num_states` is an integer and `rate` is a numeric. Creates a `num_states` x `num_states` rate matrix with average rate `rate`.
#' - `makeRateMatrix(matrix = Q)` where `Q` is a `base::matrix`. Creates a rate matrix with the same values as in `Q`.
#'
#' ## Methods
#' Call a method on a `RateMatrix` object `Q` with `Q$methodName(arguments)`.
#' - `getMatrix()`: Returns a standard R [base::matrix].
#'
#' ## Accessors
#' `RateMatrix` provides intuitive access to elements of the rate matrix, but enforces some behavior.
#' - `[i,j]` returns value in the _i_th row and _j_th column.
#' - `[i,j] <- y` sets the value of the _ij_th element to y. **_Will result in an error if you attempt to set a diagonal value, or attempt to set an off-diagonal element to a negative value._**
#'
#' @aliases RateMatrix, RateMatrixList
#'
#' @examples
#' # create a 4x4 rate matrix with average rate of 0.1 (ie a Jukes-Cantor model with mu = 0.1).
#' Q <- makeRateMatrix(4, 0.1)
#'
#' # set the rate from 1 to 2 to 0.2
#' Q[1,2] <- 0.2
#'
#' \dontrun{
#' # setting a diagonal value is prohibited
#' Q[1,1] <- 0.1
#'
#' # likewise, setting an off-diagonal value to a negative number is prohibited
#' Q[1,2] <- -0.1
#' }
#'
#' # construct a time-homogeneous rate matrix model
#' Q  <- makeRateMatrix(4, 0.1)
#' tp <- new(TensorPhyloInstance, 4)
#' tp$setEtaConstantUnequal(Q)
#'
#' # construct a time-heterogeneous rate matrix model
#' Q_1 <- makeRateMatrix(4, 0.1)
#' Q_2 <- makeRateMatrix(4, 0.2)
#' Qs  <- c(Q_1, Q_2)
#' t   <- 1
#' tp  <- new(TensorPhyloInstance, 4)
#' tp$setEtaTimeDependentUnequal(t, Qs)
#'
#' @export
makeRateMatrix <- function(num_states = NULL, rate = 0, matrix = NULL) {
  if ( is.null(matrix) == TRUE & is.null(num_states) == FALSE ) {
    Q <- new(RateMatrix, num_states, rate)
  } else if ( is.null(matrix) == FALSE ) {
    Q <- new(RateMatrix, matrix)
  } else {
    stop("Something went wrong trying to create the rate matrix. Please see the documentation for help.")
  }
  return(Q)
}
