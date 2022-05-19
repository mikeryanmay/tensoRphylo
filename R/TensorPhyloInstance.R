#' @title API documentation
#'
#' @section Constructor:
#'
#' @section Models and Parameters:
#'
#' ## Lineage-specific diversification events
#'
#' ### Speciation events
#'
#' ### Extinction events
#'
#' ### Sampling events
#'
#' ### Destructive-sampling events
#'
#' ## Mass-diversification events
#'
#' ### Mass-speciation events
#'
#' ### Mass-extinction events
#'
#' ### Mass-sampling events
#'
#' ### Mass-destructive-sampling events
#'
#' ## State-change events
#'
#' ### Anagenetic events
#'
#' ### Cladogenetic events
#'
#' ## Root Frequencies
#'
#' ## Conditional probabilities
#'
#' @section Numerical settings:
#'
#' @name TensorPhyloInstance
#' @aliases makeTensorPhylo
#' @export
NULL

#' @export
makeTensorPhylo <- function(tree, data = NULL, num_states = NULL) {

  # check the tree format
  if ( inherits(tree, "phylo") == FALSE ) {
    stop("Tree must be in ape::phylo format.")
  }

  # check that the tree is rooted
  if ( is.rooted(tree) == FALSE ) {
    stop("Tree must be rooted.")
  }

  # check that the tree is not multifurcating
  if ( .is.multifurcating(tree) == TRUE ) {
    stop("Tree must not have polytomies.")
  }

  # add a small stem to the tree, if necessary
  if ( is.null(tree$root.edge) == TRUE ) {
    tree$root.edge <- 0.0
  }

  # get the data as a newick string
  newick <- write.tree(tree)

  # if data are empty, make an empty data frame
  if ( is.null(data) == TRUE ) {

    # make sure we specified number of states
    if ( is.null(num_states) == TRUE ) {
      stop("If you don't provide data, you must provide the number of missing states.")
    }

    # make an empty data matrix
    # make sure to include sampled ancestors (labeled nodes)
    sample_labels  <- c(tree$tip.label, tree$node.label)
    sample_labels  <- sample_labels[sample_labels != ""]
    data           <- matrix(1, nrow = length(sample_labels), ncol = num_states)
    rownames(data) <- sample_labels
    colnames(data) <- 1:num_states - 1

  }

  # convert a vector to a data matrix, if necessary
  if ( is.vector(data) == TRUE ) {
    data <- .char.vector.to.table(data)
  }

  # check the data format
  if ( is.matrix(data) == FALSE ) {
    stop("Data must be provided in matrix format.")
  }

  # check consistency of tree and data
  if ( .is.consistent.data(tree, data) == FALSE ) {

    # throw a warning
    warning("Tree and data do not have the same samples in them. Adding missing data as appropriate.")

    # enforce conformity
    conformed_data <- .conform.data(tree, data)
    tree <- conformed_data$tree
    data <- conformed_data$data

    if ( .is.consistent.data(tree, data) == FALSE ) {
      stop("Something went wrong when trying to add missing data.")
    }

  }

  # create a new TP object
  tp <- new(TensorPhyloInstance, tree, newick, data)

  # return it
  return(tp)

}
