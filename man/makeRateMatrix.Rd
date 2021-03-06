% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RateMatrix.R
\name{makeRateMatrix}
\alias{makeRateMatrix}
\alias{RateMatrix,}
\alias{RateMatrixList}
\title{Rate Matrix}
\usage{
makeRateMatrix(num_states = NULL, rate = 0, matrix = NULL)
}
\description{
The \code{RateMatrix} class.
}
\details{
This class is for instantaneous rate matrices, which describe the rate of change from state \emph{i} (rows) to state \emph{j} (columns).
In tensorphylo, these matrices are used to describe "anagenetic", "asynchronous" events, i.e., changes in state that are not associated with speciation events, and occur independently among lineages.
This is directly analogous to the Q matrix of a substitution model.

Objects of the class \code{RateMatrix} have off-diagonal elements that are strictly positive and diagonal elements such that the sum of each row is zero.
Attempting to set the value of a diagonal element is prohibited, as they are computed automatically from the off-diagonal elements.
Likewise, the class enforces that off-diagonal elements must >= 0.

For models with time-dependent rate matrices, use the generic \code{c} to create a vector of rate matrices (of class \code{RateMatrixList}):
\preformatted{
Q_1 <- makeRateMatrix(3, 0.1)
Q_2 <- makeRateMatrix(3, 0.2)
Qs  <- c(Q_1, Q_2)
}
\subsection{Constructors}{

Create a new object of class \code{RateMatrix}.
\itemize{
\item \code{makeRateMatrix(num_states)} where \code{num_states} is an integer. Creates a \code{num_states} x \code{num_states} rate matrix with rate zero.
\item \code{makeRateMatrix(num_states, rate)} where \code{num_states} is an integer and \code{rate} is a numeric. Creates a \code{num_states} x \code{num_states} rate matrix with average rate \code{rate}.
\item \code{makeRateMatrix(matrix = Q)} where \code{Q} is a \code{base::matrix}. Creates a rate matrix with the same values as in \code{Q}.
}
}

\subsection{Methods}{

Call a method on a \code{RateMatrix} object \code{Q} with \code{Q$methodName(arguments)}.
\itemize{
\item \code{getMatrix()}: Returns a standard R \link[base:matrix]{base::matrix}.
}
}

\subsection{Accessors}{

\code{RateMatrix} provides intuitive access to elements of the rate matrix, but enforces some behavior.
\itemize{
\item \verb{[i,j]} returns value in the _i_th row and _j_th column.
\item \verb{[i,j] <- y} sets the value of the _ij_th element to y. \strong{\emph{Will result in an error if you attempt to set a diagonal value, or attempt to set an off-diagonal element to a negative value.}}
}
}
}
\examples{
# create a 4x4 rate matrix with average rate of 0.1 (ie a Jukes-Cantor model with mu = 0.1).
Q <- makeRateMatrix(4, 0.1)

# set the rate from 1 to 2 to 0.2
Q[1,2] <- 0.2

\dontrun{
# setting a diagonal value is prohibited
Q[1,1] <- 0.1

# likewise, setting an off-diagonal value to a negative number is prohibited
Q[1,2] <- -0.1
}

# construct a time-homogeneous rate matrix model
Q  <- makeRateMatrix(4, 0.1)
tp <- new(TensorPhyloInstance, 4)
tp$setEtaConstantUnequal(Q)

# construct a time-heterogeneous rate matrix model
Q_1 <- makeRateMatrix(4, 0.1)
Q_2 <- makeRateMatrix(4, 0.2)
Qs  <- c(Q_1, Q_2)
t   <- 1
tp  <- new(TensorPhyloInstance, 4)
tp$setEtaTimeDependentUnequal(t, Qs)

}
