\name{plotpsm}
\alias{plotpsm}
\title{Plot of posterior similarity matrix}
\usage{
plotpsm(psm.ind,
        psm.tot,
        method="complete",
        plot1.name,
        plot2.name,
        ...)
}
\arguments{
\item{psm.ind}{posterior similarity matrix of cells between 2 datasets, or within the same dataset. This is equivalent to the first item of the list of output from the function similarity_matrix().}
\item{psm.tot}{posterior similarity matrix of cells from all datasets. This is equivalent to the second item of the list of output from the function similarity_matrix().}
\item{method}{method used to reorder cells in each dataset. Default is to set this method to 'complete.}
\item{plot1.name}{name of the first output plot.}
\item{plot2.name}{name of the second output plot.}
\item{...}{Extra parameters for the plot}
}
\value{
The output contains two graphs. The first plots consider cells from all datasets when reordering. The second graph consider cells separately from different datasets; vertical and horizontal lines are used to separate cells from different datasets.
}
\description{
Posterior similarity plots. Note when plotting those similarity matrix, the plots output to the working directory folder directly. Viewing the plot in R directly has the potential to cause R session crashing in cases with large number of cells.
}
\examples{
# Plot posterior similarity matrix:

plotpsm(psm.ind = psm$psm.within,
        psm.tot = psm$psm.combined,
        method="complete",
        xlab = 'cells',
        ylab = 'cells',
        plot1.name = 'plot1',
        plot2.name = 'plot2')
}
