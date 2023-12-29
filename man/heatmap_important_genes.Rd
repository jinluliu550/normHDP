\name{heatmap_important_genes}
\alias{heatmap_important_genes}
\title{Heat-map of important genes}
\usage{
heatmap_important_genes(important_genes,
                        normHDP_post_output,
                        gene_list)
}
\arguments{
\item{important_genes}{a vector of names of important genes.}
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{gene_list}{a vector of names of genes in the data.}
}
\value{
heat-map of posterior estimated mean expressions and dispersions of important genes.
}
\description{
The function creates heat-maps of estimated mean expressions and dispersions for a list of selected genes.
}
\examples{
heatmap_important_genes(important_genes = important_genes_list,
                        normHDP_post_output = normHDP_post_MCMC,
                        gene_list = Y_genes)
}
