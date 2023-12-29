\name{cluster_summary}
\alias{cluster_summary}
\title{Summary of the number of cells in each cluster}
\usage{
cluster_summary(Z,
                data.name)
}
\arguments{
\item{Z}{the point estimate of clustering.}
\item{data.name}{a vector of names of each data.}
}
\value{
The function outputs 2 plot: plot 1 is a bar-chart to show the number of cells in each cluster. Plot 2 shows the proportion of data 1 cells and data 2 cells in each cluster, and compare with the overall proportion. Note plot 2 is only shown for cases with exactly 2 datasets.
}
\description{
The function summarizes the point estimate of clustering.
}
\examples{
# An example with 2 datasets:

cluster_summary(Z = Z_opt,
                data.name = c('data 1', 'data 2'))

}
