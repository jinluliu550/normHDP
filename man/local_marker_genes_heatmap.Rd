\name{local_marker_genes_heatmap}
\alias{local_marker_genes_heatmap}
\title{Heat-maps of the local marker genes}
\usage{
local_marker_genes_heatmap(lmg_output,
                           normHDP_post_output)
}
\arguments{
\item{lmg_output}{direct output from function local_marker_genes().}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
}
\value{
For each cluster, the function draws a heat-map to show the estimated unique parameters for the local marker genes. The heat-maps are shown on the log-scale.
}
\description{
Heat-maps of the local marker genes.
}
\examples{
local_marker_genes_heatmap(lmg_output = lmg_df,
                           normHDP_post_output = normHDP_post_MCMC)
}
