\name{latent_counts}
\alias{latent_counts}
\title{Posterior estimated latent counts}
\usage{
latent_counts(Y,
              normHDP_post_output,
              gmg_output = NULL,
              num.cores = 1,
              run.on.pc = TRUE)
}
\arguments{
\item{Y}{a list of input data.}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{gmg_output}{direct output from function global_marker_genes().}
\item{num.cores}{number of running cores.}
\item{run.on.pc}{an option to run on either the pc or server. Default is set to TRUE.}
}
\value{
\item{Y_latent_DE}{estimated latent counts with cells reordered by the clustering estimate and genes reordered by the tail probabilities of mean expressions. This is only returned if gmg_output is given.}
\item{Y_latent_DD}{estimated latent counts with cells reordered by the clustering estimate and genes reordered by the tail probabilities of dispersions. This is only returned if gmg_output is given.}
\item{matrix.output}{estimated latent counts with cells reordered by the clustering estimate. This is only returned if gmg_output is not given.}
\item{DE_number}{total number of global differentially expressed genes.}
\item{DD_number}{total number of global differentially dispersed genes.}
\item{index.solid}{index used for plotting the heat-map of estimated latent counts to separate cells from different clusters.}
\item{index.dashed}{index used for plotting the heat-map of estimated latent counts to separate cells from different datasets.}
\item{Z}{the optimized clustering.}
\item{G}{total number of genes.}
}
\description{
The function calculates the posterior estimated latent counts and reorder cells by the optimized clustering. There is an option to reorder genes as well if the global marker genes classification outcomes is provided.
}
\examples{
# An example with 10 cores. Ideally for cases with Y having high dimensions and a large number of samples, user should use as many cores as possible.

latent_counts(Y = Y,
              normHDP_post_output = normHDP_post_mcmc,
              gmg_output = gmg_df,
              num.cores = 10)
}
