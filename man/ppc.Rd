\name{ppc}
\alias{ppc}
\title{Posterior predictive check with multiple datasets}
\usage{
ppc(normHDP_output,
    Y,
    number_rep)
}
\arguments{
\item{normHDP_output}{Output from function normHDP_mcmc.}
\item{Y}{Input dataset of length D, where D is the total number of datasets.}
\item{number_rep}{Number of replicated datasets.}
}
\value{
\item{rep_Y_statistics}{Statistics of replicated datasets; a data frame with variables: index of dataset, index of gene, index of replicate, mean of log shifted counts, standard deviation of log shifted counts, dropout probabilities, and log of mean counts.}
\item{Y_statistics}{Statistics of the observed dataset; same structure as rep_Y_statistics, the index of replicate is set to zero in this case.}
}
\description{
To compute gene-wise statistics of the observed and replicated datasets
}
\examples{
# Load data
load('data/case1_mcmc.RData')
load('data/Y_linear.RData')

# Run posterior predictive check 30 times
# An example output is saved as case1_ppc_multiple in the data file
ppc(normHDP_output = case1_mcmc,
    Y = Y_linear,
    number_rep = 30)
}
