\name{global_marker_genes_heatmap}
\alias{global_marker_genes_heatmap}
\title{Heatmaps of Global Marker Genes}
\usage{
global_marker_genes_heatmap(gmg_output,
                            normHDP_post_output)
}
\arguments{
\item{gmg_output}{Output from function global_marker_genes.}
\item{normHDP_post_output}{Output from function normHDP_mcmc_post.}
}
\value{
The output contains 4 plots: The first plot show the posterior mean of mean expressions with genes reordered by tail probabilities; genes above the horizontal dashed line are the global differentially expressed genes. Plot 2 has the same structure as plot 1, but for dispersions. Plot 3 only contains global DE genes and mean expressions shown in the hea map are normalized to have zero mean across all clusters for each gene. Plot 4 has the same structure as plot 3, but for dispersions.
}
\description{
Heatmaps of global marker genes to identify potential patterns.
}
\examples{
# Load data
load('data/case1_global_marker.RData')
load('data/case1_mcmc_post.RData')

# Heatmap
global_marker_genes_heatmap(gmg_output = case1_global_marker,
                            normHDP_post_output = case1_mcmc_post)
}
