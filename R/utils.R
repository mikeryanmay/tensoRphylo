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

  # get the state labels
  split_data   <- strsplit(data, "/")
  state_labels <- sort(unique(unlist(split_data)))
  state_labels <- state_labels[state_labels %in% c("?","-") == FALSE]

  # replace ? and - with all possible states
  missing_data <- which(split_data %in% c("?","-"))
  for(i in 1:length(missing_data)) {
    split_data[[missing_data[i]]] <- state_labels
  }

  # make the matrix
  res <- matrix(0, nrow = length(data), ncol = length(state_labels))
  rownames(res) <- sample_labels
  colnames(res) <- state_labels

  # fill in 1 for consistent states
  res[cbind(rep(1:length(sample_labels), times = lengths(split_data)), match(unlist(split_data), state_labels))] <- 1

  return(res)

}

.conform.data <- function(tree, data) {

  # get tree labels
  tree_labels <- c(tree$tip.label, tree$node.label)
  tree_labels <- tree_labels[tree_labels != ""]

  # get the sample labels
  data_labels <- rownames(data)

  # make temporary object
  tree_tmp <- tree
  data_tmp <- data

  # if there are any tree labels NOT in data
  missing_from_data <- tree_labels[tree_labels %in% data_labels == FALSE]
  missing_data <- matrix(1, nrow = length(missing_from_data), ncol = ncol(data))
  rownames(missing_data) <- missing_from_data
  data_tmp <- rbind(data_tmp, missing_data)

  # if there are any data NOT in the tree
  data_tmp <- data_tmp[tree_labels,]

  return(list(tree = tree_tmp, data = data_tmp))

}










