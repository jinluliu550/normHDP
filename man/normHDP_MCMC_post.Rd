\name{normHDP_mcmc_post}
\alias{normHDP_mcmc_post}
\title{MCMC Simulation of mean expressions and dispersions}
\usage{
normHDP_mcmc_post(normHDP_output,
                  burn_in,
                  thinning,
                  number_iter,
                  Z,
                  auto.save = FALSE,
                  partial.save.name = NULL,
                  iter_update = 100,
                  Y,
                  num.cores = 4)
}
\arguments{
\item{normHDP_output}{Output from the normHDP_mcmc function}
\item{number_iter}{Total number of iterations.}
\item{thinning}{Default is to set thinning to 5; saves posterior estimates for every fifth iteration.}
\item{burn_in}{Default is to set to 3000; the first 3000 iterations are ignored.}
\item{Z}{Posterior point estimate of clustering.}
\item{Y}{Input dataset}
\item{auto.save}{An option to save the output every 100 iterations. Default is to set to FALSE.}
\item{partial.save.name}{Default is to set to NULL, an option to give a name to the output data if set auto.save = TRUE.}
\item{iter_update}{Controls for text, i.e.: if iter_update is set to 100, then the algorithm will output a message every 100 iterations.}
\item{num.cores}{Number of cores to be used by the function. Default is to use 4 cores.}
}
\value{
\item{mu_star_1_J_output}{Samples of mean expressions.}
\item{phi_star_1_J_output}{Samples of dispersion parameters.}
\item{mu_star_1_J_new}{The latest sample of mean expressions.}
\item{phi_star_1_J_new}{The latest sample of dispersions.}
\item{tilde_s_unique_new}{Sum of cross products of transformed mean expressions and dispersions.}
\item{mean_X_unique_new}{Sample mean of transformed mean expressions and dispersion.}
\item{covariance_unique_new}{Variance-covariance of transformed mean expressions and dispersions.}
\item{unique_count}{Total number of mean expressions and dispersions which are accepted.}
\item{acceptance_prob_vec}{Cumulative acceptance probabilities of simulation.}
\item{J}{Total number of occupied components.}
\item{G}{Number of genes.}
\item{D}{Number of datasets.}
\item{C}{Number of cells in each dataset.}
\item{alpha_mu_2}{Prior estimated variance of log mean expressions.}
\item{quadratic}{Prior assumption of relationship between mean expressions and dispersions on the log-scale.}
\item{b_posterior_mean}{Posterior mean of capture efficiencies.}
\item{b_posterior_mean}{Posterior mean of slope}
\item{final_iter}{Final iteration.}
\item{output_index}{Number of recorded samples in the output.}
\item{auto.save}{An option to save the output every 100 iterations. Default is to set to FALSE.}
\item{partial.save.name}{Default is to set to NULL, an option to give a name to the output data if set auto.save = TRUE.}
\item{iter_update}{Controls for text, i.e.: if iter_update is set to 100, then the algorithm will output a message every 100 iterations.}
\item{num.cores}{Number of cores to be used by the function. Default is to use 4 cores.}
\item{Z}{Posterior point estimate of clustering.}
\item{Y}{Input dataset}
\item{alpha_phi_2}{Estimated variance parameter for the slope.}
}
\description{
Due to the label switching property, we will need to run the following code to carry out posterior inferences on mean expressions and dispersions.
}
\examples{
# Load data
load('data/case1_mcmc.RData')
load('data/case1_Z_estimate.RData')
load('data/Y_linear.RData')

# Run MCMC
# An example saved as case1_mcmc_post
normHDP_mcmc_post(normHDP_output = case1_mcmc,
                  burn_in = 1000,
                  thinning = 5,
                  number_iter = 3000,
                  Z = case1_Z_estimate,
                  Y = Y_linear,
                  iter_update = 100)

}
