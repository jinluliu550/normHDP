\name{observed_counts}
\alias{observed_counts}
\title{Observed counts with reordered cells and genes}
\usage{
observed_counts(Z_CD,
                Y,
                gmg_output = NULL)
}
\arguments{
\item{Z}{Allocations, a list of length D, where D is the total number of datasets.}
\item{Y}{Observed dataset}
\item{gmg_output}{Output from function global_marker_genes. Default is to not provide this information.}
}
\value{
\item{matrix.output}{This is the output when gmg_output is not provided; matrix of observed counts with reordered cells by allocations.}
\item{index.solid}{Index of the boundaries between clusters in the matrix.}
\item{index.dashed}{Index of the boundaries between datasets in the matrix.}
\item{Y_latent_DE}{This is the output when gmg_output is provided; matrix with reordered cells and genes by tail probabilities corresponding to the mean expressions of genes.}
\item{Y_latent_DD}{This is the output when gmg_output is provided; matrix with reordered cells and genes by tail probabilities corresponding to the dispersions of genes.}
\item{DE_number}{Number of globally differentially expressed genes.}
\item{DD_number}{Number of globally differentially dispersed genes.}
}
\description{
Reorder genes and cells in the observed counts matrix.
}
\examples{
# Load data
load('data/case1_Z_estimate.RData')
load('data/Y_linear.RData')
loaD('data/case1_global_marker.RData')

# Reordered matrix of observed counts
# Output is saved as case1_observed_counts in the data file
observed_counts(Z_CD = case1_Z_estimate,
                Y = Y_linear,
                gmg_output = case1_global_marker)
}
