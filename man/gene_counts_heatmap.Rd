\name{latent_counts}
\alias{latent_counts}
\title{Estimated latent counts}
\usage{
latent_counts(mu_JG,
              phi_JG,
              Z_CD,
              beta_CD,
              Y,
              gmg_output = NULL)
}
\arguments{
\item{mu_JG}{Matrix of mean expressions.}
\item{phi_JG}{Matrix of dispersions.}
\item{Z}{Allocations, a list of length D, where D is the total number of datasets.}
\item{beta_CD}{Capture effiencies, a list of length D.}
\item{Y}{Observed dataset}
\item{gmg_output}{Output from function global_marker_genes. Default is to not provide this information.}
}
\value{
\item{matrix.output}{This is the output when gmg_output is not provided; matrix of latent counts with reordered cells by allocations.}
\item{index.solid}{Index of the boundaries between clusters in the matrix.}
\item{index.dashed}{Index of the boundaries between datasets in the matrix.}
\item{Y_latent_DE}{This is the output when gmg_output is provided; matrix with reordered cells and genes by tail probabilities corresponding to the mean expressions of genes.}
\item{Y_latent_DD}{This is the output when gmg_output is provided; matrix with reordered cells and genes by tail probabilities corresponding to the dispersions of genes.}
\item{DE_number}{Number of globally differentially expressed genes.}
\item{DD_number}{Number of globally differentially dispersed genes.}
}
\description{
Compute estimated latent counts.
}
\examples{
# Load data
load('data/case1_mcmc_post.RData')
load('data/case1_Z_estimate.RData')
load('data/Y_linear.RData')
load('data/case1_global_marker.RData')

# Chain length
length.mcmc <- length(case1_mcmc_post$mu_star_1_J_output)

# Posterior mean of mean expression with fixed allocation
case1.post.mean.mu <- Reduce("+", case1_mcmc_post$mu_star_1_J_output)/length.mcmc

# Posterior mean of dispersion with fixed allocation
case1.post.mean.phi <- Reduce("+", case1_mcmc_post$phi_star_1_J_output)/length.mcmc

# Posterior mean of capture efficiencies
case1.post.mean.beta <- case1_mcmc_post$Beta_posterior_mean

# Obtain reordered latent counts
# The output is saved as case1_estimated_latent in the data file
latent_counts(mu_JG = case1.post.mean.mu,
              phi_JG = case1.post.mean.phi,
              Z_CD = case1_Z_estimate,
              beta_CD = case1.post.mean.beta,
              Y = Y_linear,
              gmg_output = case1_global_marker)
}
