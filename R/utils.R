.is.multifurcating <- function(phy) any(table(phy$edge[,1]) > 2)

.is.consistent.data <- function(phy, data) {

  # make sure tree is phylo
  if ( inherits(phy, "phylo") == FALSE ) {
    stop("Tree must be in ape::phylo format.")
  }

  # make sure data is table
  if ( inherits(data, "matrix") == FALSE ) {
    stop("Data must be in base::matrix format.")
  }

  # check labels
  tree_labels <- c(phy$tip.label, phy$node.label)
  tree_labels <- tree_labels[tree_labels != ""]
  data_labels <- rownames(data)

  # return if sets are identical
  return(setequal(tree_labels, data_labels))

}

.char.vector.to.table <- function(data) {

  # make sure it's a vector
  if ( is.vector(data) == FALSE && is.null(names(data)) == TRUE ) {
    stop("Can only convert _named_ vectors data to base::matrix form. ")
  }

  # get sample labels
  sample_labels <- names(data)

  # make sure data are character
  mode(data) <- "character"
  split_data <- strsplit(data, "/")
  state_labels <- sort(unique(unlist(split_data)))

  # make the matrix
  res <- matrix(0, nrow = length(data), ncol = length(state_labels))
  rownames(res) <- sample_labels
  colnames(res) <- state_labels

  # fill in 1 for consistent states
  res[cbind(rep(1:length(sample_labels), times = lengths(split_data)), match(unlist(split_data), state_labels))] <- 1

  return(res)

}




