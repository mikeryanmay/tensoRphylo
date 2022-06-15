.zero.bl.to.SA.tree <- function(tree) {

  # find all the tips
  num_tips <- Ntip(tree)
  tips     <- tree$edge[,2] <= num_tips

  # figure out which are sampled ancestors
  sampled_ancestors <- tree$edge.length == 0 & tips
  sampled_tips      <- tree$edge.length >  0 & tips

  # compute the numbers
  num_sampled_ancestors <- sum(sampled_ancestors)
  num_sampled_tips      <- sum(sampled_tips)

  # make the new tree
  new_tree <- tree

  # re-index sampled tips
  old_tip_index <- tree$edge[sampled_tips, 2]
  new_tip_index <- old_tip_index
  new_tip_index[order(old_tip_index)] <- 1:num_sampled_tips
  new_tree$edge[sampled_tips, 2] <- new_tip_index

  # drop edges that are sampled ancestors
  new_tree$edge        <- new_tree$edge[sampled_ancestors == FALSE,]
  new_tree$edge.length <- new_tree$edge.length[sampled_ancestors == FALSE]

  # make the tip labels
  new_tree$tip.label <- new_tree$tip.label[1:num_tips %in% tree$edge[sampled_tips,2]]

  # make node labels
  new_tree$node.label <- character(new_tree$Nnode)
  new_tree$node.label[tree$edge[sampled_ancestors,1] - tree$Nnode - 1] <- tree$tip.label[tree$edge[sampled_ancestors,2]]

  # re-index nodes
  new_tree$edge[new_tree$edge > num_tips] <- new_tree$edge[new_tree$edge > num_tips] - num_sampled_ancestors

  # done
  return(new_tree)

}

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










