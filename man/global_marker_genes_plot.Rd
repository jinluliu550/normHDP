\name{global_marker_genes_plot}
\alias{global_marker_genes_plot}
\title{Summary of Global Marker Genes}
\usage{
global_marker_genes_plot(gmg_output)
}
\arguments{
\item{gmg_output}{Output from function global_marker_genes.}
}
\value{
Returns three plots: The first plot shows the relationship between absolute log-fold change and tail probabilities corresponding to mean expressions. Genes above the horizontal line are classified as global differentially expressed genes, and vice versa. Plot 2 has the same sturcture as plot 1, but for dispersions. Plot 3 summarizes the number of genes with single and both global characteristics.
}
\description{
Summary of global marker genes.
}
\examples{
# Load data
load('data/case1_global_marker.RData')

# Plot
global_marker_genes_plot(gmg_output = case1_global_marker)
}
