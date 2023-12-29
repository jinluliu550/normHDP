\name{consensus_clustering}
\alias{consensus_clustering}
\title{Finding the optimized clustering and the posterior similarity matrix.}
\usage{
consensus_clustering(Z_output,
                     num.cores,
                     run.on.pc = TRUE)
}
\arguments{
\item{Z_output}{posterior samples of allocations.}
\item{num.cores}{total number of running cores.}
\item{run.on.pc}{an option to run the algorithm on pc. Default is set to TRUE.}
}
\value{
\item{psm}{a list of length 2. The first item is labelled as 'psm.within' which is the posterior similarity matrix between 2 datasets or within the same dataset. For example: psm.within\code{[[i]][[j]]} is the similarity between cells in dataset \code{i} and \code{j}. Note that to save memory space, \code{i} is never greater than \code{j} to avoid repeats. The second item is labelled as 'psm.combined' which measure the posterior similarity of cells from all datasets.}
\item{cluster.estimate}{point estimate of clustering which minimizes the variation of information.}
}
\description{
The function calculates the posterior similarity matrix and the point estimate of the allocations.
}
\examples{
# An example with 10 running cores:

consensus_clustering(Z_output = Z_trace,
                     num.cores = 10)
}
