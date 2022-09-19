\name{local_marker_genes_plot}
\alias{local_marker_genes_plot}
\title{Plots of local marker genes}
\usage{
local_marker_genes_plot(lmg_output)
}
\arguments{
\item{lmg_output}{Output from function local_marker_genes.}
}
\value{
Output contains 4 plots. Plot 1: relationship between absolute log-fold change and tail probabilities for mean expressions with threshold shown. Plot 2: A bar chart to show number of local differentially expressed genes for each cluster. Plot 3 and 4 are similar to the structure of the previous 2 plots, but to show relationship and counts of locally differentially dispersed genes.
}
\description{
Plots of local marker genes.
}
\examples{
# Load data
load('data/case1_local_marker.RData')

# Plot
local_marker_genes_plot(lmg_output = case1_local_marker)
}
