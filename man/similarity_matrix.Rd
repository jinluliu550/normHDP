\name{similarity_matrix}
\alias{similarity_matrix}
\title{Posterior similarity matrix}
\usage{
similarity_matrix(normHDP_output)
}
\arguments{
\item{normHDP_output}{Output from function normHDP_mcmc.}
}
\value{
\item{psm.within}{posterior similarity matrix between 2 datasets or within the same dataset. For example: psm.within\code{[[i]][[j]]} is the similarity between cells in dataset \code{i} and \code{j}. Note that to save memory space, \code{i} is never greater than \code{j} to avoid repeats.}
\item{psm.combined}{Posterior similarity matrix of all cells.}
}
\description{
Compute posterior similarity matrix using posterior samples of allocations.
}
\examples{
# Load data
load('data/case1_psm.RData')

# Obtain posterior similarity matrix
similarity_matrix(normHDP_output = case1_psm)
}
