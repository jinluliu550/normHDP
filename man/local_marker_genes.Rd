\name{local_marker_genes}
\alias{local_marker_genes}
\title{Classification of local marker genes}
\usage{
local_marker_genes(normHDP_post_output,
                   threshold,
                   alpha_M = NULL,
                   alpha_D = NULL)
}
\arguments{
\item{normHDP_post_output}{Output from function normHDP_mcmc_post.}
\item{threshold}{Threshold for the absolute log-fold change between mean expressions or dispersions.}
\item{alpha_M}{Genes with tail probabilities of mean expressions greater than alpha_M are classified as local differentially expressed genes. By default, this is set to control expected false discovery rate to 5 percent.}
\item{alpha_D}{Genes with tail probabilities of dispersions greater than alpha_D are classified as local differentially dispersed genes. By default, this is set to control expected false discovery rate to 5 percent.}
}
\value{
\item{marker_DE}{A data frame with 5 variables; including cluster, gene, tail probability, absolute log-fold change and a binary variable to indicate whether the gene is a local differentially expressed gene.}
\item{marker_DD}{A data frame with 5 variables; including cluster, gene, tail probability, absolute log-fold change and a binary variable to indicate whether the gene is a local differentially dispersed gene.}
\item{alpha_M}{Threshold for the tail probability of mean expressions.}
\item{alpha_D}{Threshold for the tail probability of dispersions.}
}
\description{
Classification of local marker genes with corressponding absolute log-fold change and tail probabilities.
}
\examples{
# Load data
load('data/case1_mcmc_post.RData')

# Obtain local marker genes
# An example is saved as case1_local_marker in the data file
local_marker_genes(normHDP_post_output = case1_mcmc_post,
                                         threshold = list('mu' = 1, 'phi' = 1),
                                         alpha_M = 0.9,
                                         alpha_D = 0.9)
}
