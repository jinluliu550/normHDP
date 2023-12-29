\name{beta_compare}
\alias{beta_compare}
\title{Comaprison of capture efficiencies}
\usage{
beta_compare(Y,
             beta.mean = 0.06,
             normHDP_post_output,
             data.names)
}
\arguments{
\item{Y}{the input dataset.}
\item{beta.mean}{prior mean of capture efficiencies, default is set to 0.06.}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixd_z().}
\item{data.names}{vector of names of each data.}
}
\value{
The function output 2 graphs: the first graph compares the baynorm estimated capture effiencices with the posterior means using scatter plot, for each dataset. The second graph compares the overall distribution of capture effiencices using histograms, for each dataset.
}
\description{
Comparison of the baynorm estimated capture efficiencies with the posterior means.
}
\examples{
# An example with 2 datasets:

beta_compare(Y = Y,
             normHDP_post_output = normHDP_post_MCMC,
             data.names = c('data 1', 'data 2'))
}
