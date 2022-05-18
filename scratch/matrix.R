library(expm)


Q <- makeRateMatrix(3)
Q <- makeRateMatrix(3, 0.1)

B <- matrix(2, nrow = 2, ncol = 2)
diag(B) <- -2

Q <- new(RateMatrix, B)
