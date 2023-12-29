\name{ppc_single_plot}
\alias{ppc_single_plot}
\title{(mixed) Posterior predictive checks with single set of replicated data}
\usage{
ppc_single_plot(normHDP_post_output,
                Y,
                data.name)
}
\arguments{
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{Y}{a list of input data.}
\item{data.name}{a vector of name of each data}
}
\value{
The output consists of 2 parts: part 1 compares 1) the relationship between mean and standard deviation of log shifted counts, and 2) the relationship beween log of mean counts and dropout probabilities, for the replicated and the observed counts of each data. Part 2) compare three statistics between the observed and replicated dataset which includes the mean and standard deviation of log shifted counts and the dropout probabilities.
}
\description{
The function compare the observed data with a single set of replicated data.
}
\examples{
# An example with 2 datasets:

ppc_single_plot(normHDP_post_output = normHDP_post_MCMC,
                Y = Y,
                data.name = c('data 1', 'data 2'))
}
