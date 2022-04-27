library(expm)


M <- new(RateMatrix, 4, 0.1)
M$setRate(1,2, 3.0)

Q <- M$getMatrix()
