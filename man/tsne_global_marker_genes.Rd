\name{tsne_global_marker_genes}
\alias{tsne_global_marker_genes}
\title{t-SNE plot with global marker genes}
\usage{
tsne_global_marker_genes(gmg_output,
                         Y,
                         Y_latent,
                         normHDP_post_output)
}
\arguments{
\item{Y}{a list of input data.}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{gmg_output}{direct output from function global_marker_genes().}
\item{Y_latent}{posterior estimated latent counts}
}
\value{
The function outputs 2 plots. Plot 1 is the t-SNE plot of the posterior estimated latent counts with global marker genes. Plot 2 is the t-SNE plot of the observed counts with global marker genes. Note only genes which are classified as being both global DE and DD are selected.
}
\description{
The function creates t-SNE plot of the observed and posterior estimated counts with global marker genes.
}
\examples{
tsne_global_marker_genes(gmg_output = gmg_df,
                         Y = Y,
                         Y_latent = Y_latent,
                         normHDP_post_output = normHDP_MCMC)
}
