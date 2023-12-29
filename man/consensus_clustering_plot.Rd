\name{consensus_clustering_plot}
\alias{consensus_clustering_plot}
\title{Effect of width and length on the mean absolute difference in posterior similarity matrix}
\usage{
consensus_clustering_plot(Z_output,
                          W.breaks,
                          L.breaks,
                          num.cores,
                          run.on.pc)
}
\arguments{
\item{Z_output}{posterior samples of allocations.}
\item{W.breaks}{a vector of breaks of chain width.}
\item{L.breaks}{a vector of breaks of chain length.}
\item{num.cores}{total number of running cores.}
\item{run.on.pc}{an option to set either or not to run the function on PC. Default is set to TRUE.}
}
\value{
A line graph to show the effect of changing chain widths and lengths on the mean absolute difference in posterior similarity matrix.
}
\description{
The function calculates the difference in posterior similarity matrix at different chain widths and depths.
}
\examples{
# An example with 10 running cores:

consensus_clustering_plot(Z_output = Z_trace,
                          W.breaks = c(10,30,50,70,90,100),
                          L.breaks = c(10,20,30,40,50,60,70,80,90,100),
                          num.cores = 10)
}
