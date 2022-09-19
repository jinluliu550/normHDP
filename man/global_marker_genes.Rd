\name{global_marker_genes}
\alias{global_marker_genes}
\title{Classification of Global Marker Genes}
\usage{
global_marker_genes(normHDP_post_output,
                    threshold,
                    alpha_M = NULL,
                    alpha_D = NULL)
}
\arguments{
\item{normHDP_post_output}{Output from function normHDP_mcmc_post.}
\item{threshold}{Threshold value for absolute log-fold change. Values are given for both mean expression and dispersion.}
\item{alpha_M}{Threshold for the tail probabilities corresponding to the mean expressions, genes with tail probabilities greater than alpha_M are classified as global differentially expressed genes. By default, this is set to control the expected false discovery rate to 5 percent.}
\item{alpha_D}{Threshold for the tail probabilities corresponding to dispersions, genes with tail probabilities greater than alpha_D are classified as global differentially dispersed genes. By default, this is set to control the expected false discovery rate to 5 percent.}
}
\value{
\item{marker_DE}{A data frame with 4 variables; including gene index, absolute log-fold change, tail-probabilities and an indicator to indicate whether the gene is a global differentially expressed gene.}
\item{marker_DD}{A data frame with 4 variables; including gene index, absolute log-fold change, tail-probabilities and an indicator to indicate whether the gene is a global differentially dispersed gene.}
\item{alpha_M}{A threshold for the tail probability corresponding to the mean expressions.}
\item{alpha_D}{A threshold for the tail probability corresponding to the dispersions.}
}
\description{
Classification of global marker genes.
}
\examples{
# Load data
load('data/case1_mcmc_post.RData')

# Obtain global marker genes based on the threshold below
# Result is saved as case1_global_marker in the data file
global_marker_genes(normHDP_post_output = case1_mcmc_post,
                    threshold = list('mu' = 1.5, 'phi' = 1.5),
                    alpha_M = 0.9,
                    alpha_D = 0.9)
}
