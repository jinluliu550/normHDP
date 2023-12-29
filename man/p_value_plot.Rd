\name{p_value_plot}
\alias{p_value_plot}
\title{Distribution of p-values of posterior predictive checks}
\usage{
p_value_plot(ppc_output)
}
\arguments{
\item{ppc_output}{direct output from function ppp_calculation().}
}
\value{
histogram of p-values for each of the three discrepancy measures, namely, 1) the median difference between the observed and expected Y for each gene, 2) the IQR of the difference between the observed and expected Y for each gene, and 3) expected difference between the true and expected dropout probabilities for each gene.
}
\description{
The function draws histograms of p-values of posterior predictive checks.
}
\examples{
p_value_plot(ppc_output = ppc_df)
}
