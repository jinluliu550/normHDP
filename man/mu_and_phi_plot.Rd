\name{mu_and_phi_plot}
\alias{mu_and_phi_plot}
\title{Posterior relationship between mean expressions and dispersions}
\usage{
mu_and_phi_plot(normHDP_post_output,
                row.n)
}
\arguments{
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{row.n}{number of rows for the plot.}
}
\value{
For each cluster, the function generates a plot to show relationship between estimated mean expressions and dispersions on the log-scale. Each graph also shows the 95 percent credible bond of the posterior relationship.
}
\description{
The function to show the posterior relationship between mean expressions and dispersions.
}
\examples{
# An example with 4 rows:

mu_and_phi_plot(normHDP_post_output = normHDP_MCMC_post,
                row.n = 4)
}
