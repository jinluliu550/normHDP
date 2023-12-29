\name{tsne_plot}
\alias{tsne_plot}
\title{t-SNE plot of observed and estimated counts (using all available genes)}
\usage{
tsne_plot(normHDP_post_output,
          Y_latent,
          Y)
}
\arguments{
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{Y_latent}{estimated latent counts.}
\item{Y}{observed counts.}
}
\value{
The function returns 2 scatter plots: Plot 1 is the t-SNE plot of the estimated latent counts. Plot 2 is the t-SNE plot for the observed counts.
}
\description{
t-SNE plot for the observed counts and estimated latent counts, colour-coded by the optimized clustering.
}
\examples{
tsne_plot(normHDP_post_output = normHDP_post_MCMC,
          Y_latent = Y_latent_df,
          Y = Y)
}
