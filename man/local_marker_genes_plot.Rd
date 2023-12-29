\name{local_marker_genes_plot}
\alias{local_marker_genes_plot}
\title{Graph for the local marker genes}
\description{
The function outputs graphs for the local marker genes
}
\usage{
local_marker_genes_plot(lmg_output)
}
\arguments{
\item{lmg_output}{direct output from function local_marker_genes().}
}
\value{
The function outputs four plots: plot 1 and 3 shows the scatter plot of absolute log-fold change against the tail probabilities corresponding to the mean expressions and dispersions, respectively, for each cluster. Local marker genes are shown in red and non-marker genes are shown in grey. Plot 2 and 4 are the bar-charts to summarize the total number of local marker genes for each cluster.
}
\examples{
local_marker_genes_plot(lmg_output = lmg_data)
}
