\name{global_marker_genes_plot}
\alias{global_marker_genes_plot}
\title{Graph for the global marker genes}
\description{
The function outputs graphs for the global marker genes
}
\usage{
global_marker_genes_plot(gmg_output)
}
\arguments{
\item{gmg_output}{direct output from function global_marker_genes().}
}
\value{
The function outputs three graphs: graph 1 and 2 plots the aboslute log-fold change against tail probabilities corressponding to the mean expressions and dispersions, respectively. Graph 3 is a bar-chart summarizing the number of each type of global marker genes.
}
\examples{
global_marker_genes_plot(gmg_output = gmg_data)
}
