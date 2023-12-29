\name{global_marker_genes_heatmap}
\alias{global_marker_genes_heatmap}
\title{Heat-maps of the global marker genes}
\usage{
global_marker_genes_heatmap(gmg_output,
                            normHDP_post_output)
}
\arguments{
\item{gmg_output}{direct output from function global_marker_genes().}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
}
\value{
The function outputs 4 graphs: graph 1 and 2 are the heat-maps for the estimated mean expressions and dispersions, respectively, with genes reordered by the tail probabilities. Global marker and non-global marker genes are separated by a horizontal line. Graph 3 and 4 are the heat-maps with only the global marker genes, genes are also reordered by the tail probabilities and standardized to have the mean of unique parameters for each gene equal to zero.
}
\description{
Heat-maps of the global marker genes.
}
\examples{
global_marker_genes_heatmap(gmg_output = gmg_df,
                            normHDP_post_output = normHDP_post_MCMC)
}
