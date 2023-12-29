\name{opt.clust}
\alias{opt.clust}
\title{Compute the cluster estimate using all allocation samples}
\usage{
opt.clust(Z_trace,
          psm_output)
}
\arguments{
\item{Z_trace}{posterior samples of allocations.}
\item{psm_output}{direct output from function similarity_matrix().}
}
\value{
Calculate the cluster point estimate which minimizes the variation of information.
}
\description{
The algorithm considers all allocation samples and the posterior similarity matrix to compute a point estimate which minimizes the variation of information.
}
\examples{

opt.clust(Z_trace = Z_sample,
          psm_output = psm)

}
