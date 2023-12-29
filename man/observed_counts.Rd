\name{observed_counts}
\alias{observed_counts}
\title{Reordering of observed counts}
\usage{
observed_counts(normHDP_post_output,
                Y,
                gmg_output = NULL)
}
\arguments{
\item{Y}{a list of input data.}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{gmg_output}{direct output from function global_marker_genes().}
}
\value{
\item{Y_DE}{observed counts with cells reordered by the clustering estimate and genes reordered by the tail probabilities of mean expressions. This is only returned if gmg_output is given.}
\item{Y_DD}{observed counts with cells reordered by the clustering estimate and genes reordered by the tail probabilities of dispersions. This is only returned if gmg_output is given.}
\item{matrix.output}{observed counts with cells reordered by the clustering estimate. This is only returned if gmg_output is not given.}
\item{DE_number}{total number of global differentially expressed genes.}
\item{DD_number}{total number of global differentially dispersed genes.}
\item{index.solid}{index used for plotting the heat-map of observed counts to separate cells from different clusters.}
\item{index.dashed}{index used for plotting the heat-map of observed counts to separate cells from different datasets.}
\item{Z}{the optimized clustering.}
}
\description{
The function reorders observed counts by the optimized clustering and also has an additional function to reorder counts by tail probabilities.
}
\examples{
observed_counts(Y = Y,
                normHDP_post_output = normHDP_post_mcmc,
                gmg_output = gmg_df)
}
