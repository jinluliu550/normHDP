\name{ppc_single}
\alias{ppc_single}
\title{Posterior predictive check with single replicate}
\usage{
ppc_single(normHDP_output,
           Y,
           data.name)
}
\arguments{
\item{normHDP_output}{Output from function normHDP_mcmc.}
\item{Y}{Input dataset; a list of length D, where D is the total number of datasets.}
\item{data.name}{Name of each dataset, a vector of length D.}
}
\value{
For each dataset, the function outputs plot to compare relationship between mean of log shifted counts and standard deviation of log shifted counts for the observed dataset against the replicated dataset, and also compare relationship between log of mean counts and dropout probabilities. In addition, statistics are compared gene-wise, with x-axis showing the statistics from the observed data and y-axis showing the statistics from the replicated data. If the model fits the data well, then these statistics between the replicated and the observed should be similar.
}
\description{
To carry out posterior predictive checks with single replicate of data.
}
\examples{
# Load data
load('data/case1_mcmc.RData')

# Repeat the single posterior predictive check twice
for(i in 1:2){

  ppc_single(normHDP_output = case1_mcmc,
             Y = Y_linear,
             data.name = c('data 1', 'data 2'))
}

}
