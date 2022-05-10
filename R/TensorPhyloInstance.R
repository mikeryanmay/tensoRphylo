#' Interface
#' @description The tensorphylo interface.
#'
#' @details
#' This is the general interface for working with the tensorphylo high-performance library for state-dependent birth-death models.
#'
#' @section Creating TensorPhyloInstance objects:
#'
#' ## Hmm
#'
#' @section Speciation rates
#' @section Extinction rates
#' @section Sampling rates
#' @section Destructive-sampling rates
#' @section Mass-speciation events
#' @section Mass-extinction events
#' @section Mass-sampling events
#' @section Mass-destructive-sampling events
#' @section Anagenetic state-change events
#' @section Cladogenetic state-change events
#' @section Mass-anagenetic-state-change events
#'
#' @name Interface
#' @aliases TensorPhylo
#' @export
makeTensorPhylo <- function(tree, data = NULL, nstates = NULL) {

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
    if ( is.null(nstates) == TRUE ) {
      stop("If you don't provide data, you must provide the number of missing states.")
    }

    # make an empty data matrix
    # make sure to include sampled ancestors (labeled nodes)
    sample_labels  <- c(tree$tip.label, tree$node.label)
    sample_labels  <- sample_labels[sample_labels != ""]
    data           <- matrix(1, nrow = length(sample_labels), ncol = nstates)
    rownames(data) <- sample_labels
    colnames(data) <- 1:nstates - 1

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
    stop("Tree and data do not have the same samples in them.")
  }

  # create a new TP object
  tp <- new(TensorPhyloInstance, tree, newick, data)

  # return it
  return(tp)

}
