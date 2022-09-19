\name{opt.clust}
\alias{opt.clust}
\title{Point estimate of allocations}
\usage{
opt.clust(normHDP_output,
          psm_output)
}
\arguments{
\item{normHDP_output}{Output from function normHDP_mcmc.}
\item{psm_output}{Output from function similarity_matrix.}
}
\value{
Point estimate of allocations, a list of length D, where D is the total number of datasets.
}
\description{
To obtain a point estimate of optimal clustering which minimizes the Variation of Information using all allocation samples.
}
\examples{
# Load data
load('data/case1_mcmc.RData')
load('data/case1_psm.RData')

# Obtain posterior point estimate of Z
# Saved as case1_Z_estimate in the data file
opt.clust(normHDP_output = case1_mcmc,
          psm_output =  case1_psm)
}
