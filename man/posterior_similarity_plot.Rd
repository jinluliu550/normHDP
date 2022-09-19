\name{plotpsm}
\alias{plotpsm}
\title{Plot of posterior similarity matrix}
\usage{
plotpsm(psm.ind,
        psm.tot,
        method="complete",
        ...)
}
\arguments{
\item{psm.ind}{Posterior similarity matrix of cells between 2 datasets, or within the same dataset. This is equivalent to the first item of the output from the function similarity_matrix.}
\item{psm.tot}{Posterior similarity matrix of cells between all datasets. This is equivalent to the second item of the output from the function similarity_matrix.}
\item{method}{Method used to reorder cells in each dataset. Default is to set this method to 'complete.}
\item{...}{Extra parameters for the plot}
}
\value{
The output contains two graphs. The first graph consider cells from all datasets when reordering. The second graph consider cells separately from different datasets; vertical and horizontal lines are used to separate cells from different datasets.
}
\description{
Posterior similarity plots.
}
\examples{
# Load data
load('data/case1_psm.RData')

# Plot posterior similarity matrix
plotpsm(psm.ind = case1_psm$psm.within,
        psm.tot = case1_psm$psm.combined,
        method="complete",
        xlab = 'cells',
        ylab = 'cells')
}
