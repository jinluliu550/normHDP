\name{ppc_multiple}
\alias{ppc_multiple}
\title{(mixed) Posterior predictive checks with multiple replicated datasets}
\usage{
ppc_multiple(normHDP_post_output,
             Y,
             number_rep,
             num.cores,
             run.on.pc = TRUE)
}
\arguments{
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{Y}{a list of input data.}
\item{number_rep}{number of replicated datasets.}
\item{num.cores}{number of running cores.}
\item{run.on.pc}{an option to run the algorithm on a PC. Default is set to TRUE.}
}
\value{
\item{rep_Y_statistics}{statistics of replicated datasets.}
\item{Y_statistics}{statistics of the observed dataset.}
\item{number_rep}{number of replicated datasets.}
}
\description{
The function calculates the statistics for the observed data and each of the replicated datasets.
}
\examples{
# An example with 20 replicates running with 10 cores:

ppc_multiple(normHDP_post_output = normHDP_post_MCMC,
             Y = Y,
             number_rep = 20,
             num.cores = 10)
}
