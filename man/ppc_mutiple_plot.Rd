\name{ppc_multiple_plot}
\alias{ppc_multiple_plot}
\title{Plot for the (mixed) posterior predictive checks with multiple replicated datasets.}
\usage{
ppc_multiple_plot(ppc_multiple_df,
                  title)
}
\arguments{
\item{ppc_multiple_df}{direct output from function ppc_multiple().}
\item{title}{a vector of names of each dataset.}
}
\value{
Kernel density plots to compare the distribution of statistics between the observed and each of the replicated datasets.
}
\description{
The function compares the kernel density estimates of mean and standard deviation of log shifted counts and dropout probabilities between the observed and each replicated dataset. The kernel density estimate of the observed data is shown in red and the ones for the replicated data are shown in grey.
}
\examples{
# An example with 2 sets of data:

ppc_multiple_plot(ppc_multiple_df = ppc_multiple_output,
                  title = c('data 1', 'data 2'))
}
