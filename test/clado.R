library(abind)
library(Matrix)
library(slam)
library(tensorr)

# just a simple test
tp <- new(TensorPhylo, 2)
tp$setDebugMode(1)

# make a sparse matrix
dim <- 2
O <- array(dim = c(dim,dim,dim))

O_sparse <- as.simple_sparse_array(O)
O_sparse[1,1,1] <- 0.5

reduce_simple_sparse_array(O_sparse)

# O[1,1,1] <- 0.50
# O[1,1,2] <- 0.25
# O[1,2,1] <- 0.25
# O[2,2,2] <- 0.50
# O[2,1,2] <- 0.25
# O[2,2,1] <- 0.25



O <- sptensor(matrix(c(NA,NA,NA), ncol= 1), vals = NA, dims = c(2,2,2))
O[1,1,1] <- 0.5
O[1,1,2] <- 0.25
O[1,2,1] <- 0.25

O@subs
O@vals
