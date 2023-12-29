\name{beta_analysis}
\alias{beta_analysis}
\title{Summary of capture efficiencies}
\usage{
beta_analysis(normHDP_post_output,
              cluster.labels)
}
\arguments{
\item{normHDP_post_output}{output from function normHDP_MCMC_fixed_z().}
\item{cluster.labels}{output from function difference_in_p_jd2().}
}
\value{
Point estimate of allocations, a list of length D, where D is the total number of datasets.
}
\description{
Draw box-plots to summarize capture efficiencies of cells in each cluster, colour-coded by the type of cluster.
The function should only be used in cases with 2 datasets.
}
\examples{

beta_analysis(normHDP_post_output = normHDP_post_MCMC,
              cluster.labels = difference_in_pjd2$df0_all)

}
