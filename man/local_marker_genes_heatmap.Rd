\name{local_marker_genes_heatmap}
\alias{local_marker_genes_heatmap}
\title{Heatmap of local marker genes}
\usage{
local_marker_genes_heatmap(lmg_output,
                           normHDP_post_output)
}
\arguments{
\item{lmg_output}{Output from function local_marker_genes.}
\item{normHDP_post_output}{Output from function normHDP_mcmc_post.}
}
\value{
Output contains 2 times the total number of clusters plots. For mean expressions, for each cluster, show a heat-map for locally differentially expressed genes. Genes are reordered by ascending tail probabilities. Graphs for dispersions have the same structure as for mean expressions.
}
\description{
Heatmap of local marker genes.
}
\examples{
# Load data
load('data/case1_local_marker.RData')
load('data/case1_mcmc_post.RData')

# Heatmap
local_marker_genes_heatmap(lmg_output = case1_local_marker,
                           normHDP_post_output = case1_mcmc_post)
}
