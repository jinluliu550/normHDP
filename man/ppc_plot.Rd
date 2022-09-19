\name{ppc_plot}
\alias{ppc_plot}
\title{Plots for comparing statistics between replicated and the observed dataset}
\usage{
ppc_plot(ppc_output,
         data.name)
}
\arguments{
\item{ppc_output}{Output from function ppc.}
\item{data.name}{Name of each dataset; a vector of length D, where D is the total number of datasets.}
}
\value{
For each dataset, plot the fitted relationship between mean of log shifted counts against standard deviation of log shifted counts, and plot the fitted relationship between log of mean shifted counts against dropout probabilities. Fitted relationship for the observed data are plotted in red. Fitted relationship for the replicated data are plotted in grey.
}
\description{
Plots to compare statistics of replicated datasets and the observed dataaset.
}
\examples{
# Load data
load('data/case1_mcmc_multiple.RData')

# Run multiple posterior predictive checks
ppc_plot(ppc_output = case1_ppc_multiple,
         data.name = c('data1', 'data2'))
}
