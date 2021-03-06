% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readData.R
\name{readDelimitedData}
\alias{readDelimitedData}
\title{Read character data in a delimited text}
\arguments{
\item{file}{path to a delimited text file (see details for format).}

\item{delim}{the delimiter, either ",", or "\\t".}

\item{nstates}{the number of states.}
}
\value{
A matrix formatted for input to tensorphylo. Named rows correspond to species, columns correspond to character states. A value of 1 indicates that the state is consistent with the observed data.
}
\description{
Read character data in a delimited text
}
\details{
Read data in a delimited text file.
Useful for cases where nexus doesn't accommodate your dataset (e.g., discrete characters with lots of states, like chromosome numbers).
The file may be comma- or tab-delimited, with the first column corresponding to the species name, and the second column corresponding to the character state(s).
States should be encoded as numbers between 0 and \code{nstates - 1} (i.e., the first state is 0).
Ambiguous data can be encoded by listing all consistent states, separated by \code{/}. For example, \code{2/3/7} indicates the species is possibly in state 2, 3 or 7.
Question marks, \verb{?}, or dashes ,\code{-}, indicate missing data.
}
\examples{
# path to character data file
data_file <- system.file("testdata", "sampled_ancestor_data.csv",
                         package = "tensoRphylo")

# read the character data
data <- readDelimitedData(data_file, delim = ",", nstates = 2)

data
}
