#' @export
tensoRphylo <- function(tree, data = NULL, nstates = NULL) {

  # add a small stem to the tree, if necessary
  if ( is.null(tree$root.edge) ) {
    tree$root.edge <- 0.0
  }

  # get the data as a newick string
  newick <- write.tree(tree)

  # if data are empty, make an empty data frame
  if ( is.null(data) ) {

    # make sure we specified number of states
    if ( is.null(nstates) ) {
      stop("If you don't provide data, you must provide the number of missing states.")
    }

    # make an empty data matrix
    sample_labels <- c(tree$tip.label, tree$node.label)
    data <- matrix(1, nrow = length(sample_labels), ncol = nstates)
    rownames(data) <- sample_labels
    colnames(data) <- 1:nstates - 1
    # colnames(data) <- 1:nstates

  }

  # create a new TP object
  tp <- new(TensorPhylo, tree, newick, data)

  # return it
  return(tp)

}
