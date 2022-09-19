\name{normHDP_mcmc_2}
\alias{normHDP_mcmc_2}
\title{MCMC Simulation - Continuation}
\usage{
normHDP_mcmc_2(normHDP_output,
               number_iter)
}
\arguments{
\item{normHDP_output}{Output from the normHDP_mcmc function.}
\item{number_iter}{Total number of iterations to run.}
}
\value{
\item{b_output}{Samples of estimated slope parameters.}
\item{alpha_phi_2_output}{Samples of estimated variance parameters for the slope.}
\item{Z_output}{Samples of allocations.}
\item{Z_new}{The latest allocation sample.}
\item{P_J_D_output}{Samples of dataset specific component probabilities.}
\item{P_J_D_new}{The latest sample of dataset specific component probabilities.}
\item{P_output}{Samples of component probabilities.}
\item{P_new}{The lastest sample of component probabilites.}
\item{P_count}{Total number of component probabilities which are accepted.}
\item{tilde_s_component_new}{Sum of cross products of the transformed component probabilities.}
\item{mean_X_component_new}{Sample mean of the transformed component probabilities.}
\item{covariance_component_new}{Sample covariance of the transformed component probabilities.}
\item{allocation_prob_output}{Samples of allocation probabilities of cells to the allocated clusters in each sample.}
\item{alpha_output}{Samples of alpha which controls for variation of dataset specific component probabilities between datasets.}
\item{alpha_new}{The latest sample of alpha.}
\item{mean_X_alpha_new}{Sample mean of transformed alpha.}
\item{M_2_alpha_new}{Sum of square differences between observation and mean of transformed alpha.}
\item{alpha_count}{Total number of alpha which are accepted.}
\item{variance_alpha_new}{Sample variance of transformed alpha.}
\item{alpha_zero_output}{Samples of alpha_zero which controls for variation of component probabilities.}
\item{alpha_zero_new}{The latest sample of alpha_zero.}
\item{mean_X_alpha_zero_new}{Sample mean of transformed alpha_zero.}
\item{M_2_alpha_zero_new}{Sum of square differences between observation and mean of transformed alpha_zero.}
\item{variance_alpha_zero_new}{Sample variance of transformed alpha_zero.}
\item{alpha_zero_count}{Total number of alpha_zero which are accepted.}
\item{mu_star_1_J_output}{Samples of mean expressions.}
\item{mu_star_1_J_new}{The latest sample of mean expressions.}
\item{phi_star_1_J_output}{Samples of dispersion parameters.}
\item{phi_star_1_J_new}{The latest sample of dispersions.}
\item{tilde_s_unique_new}{Sum of cross products of transformed mean expressions and dispersions.}
\item{mean_X_unique_new}{Sample mean of transformed mean expressions and dispersion.}
\item{covariance_unique_new}{Variance-covariance of transformed mean expressions and dispersions.}
\item{unique_count}{Total number of mean expressions and dispersions which are accepted.}
\item{Beta_output}{Samples of capture efficiencies (Beta).}
\item{Beta_new}{The lastest sample of Beta.}
\item{mean_X_capture_new}{Sample mean of transformed Beta.}
\item{M_2_capture_new}{Sum of square differences between observation and sample mean of Beta.}
\item{variance_capture_new}{Sample variance of Beta.}
\item{Beta_count}{Total counts of accepted Beta.}
\item{output_index}{Total number of samples in the output.}
\item{acceptance_prob_list}{Data frame to capture acceptance probabilities of parameters.}
\item{alpha_mu_2}{Prior estimated variance of log mean expressions.}
\item{a_d_beta}{Prior estimated shape parameter 1 of Beta.}
\item{b_d_beta}{Prior estimated shape parameter 2 of Beta.}
\item{quadratic}{Prior assumption of relationship between mean expressions and dispersions on the log-scale.}
\item{iter}{Final iteration number.}
\item{adaptive_prop}{Additional variance of the random walk.}
\item{J}{Total number of components.}
\item{G}{Number of genes.}
\item{D}{Number of datasets.}
\item{C}{Number of cells in each dataset.}
\item{thinning}{saves posterior estimates for every specified iteration.}
\item{empirical}{Method to estimated fixed parameters in priors.}
\item{burn_in}{First specified number of iterations to ignore.}
\item{iter_update}{Controls for text, i.e.: if iter_update is set to 100, then the algorithm will output a message every 100 iterations.}
\item{partial.save.name}{Default is to set to NULL, an option to give a name to the output data if set auto.save = TRUE.}
\item{adpative_prop}{Extra variation of the random walk.}
\item{BB_SIZE}{An option to give to the bayNorm function to determine the method used to estimate dispersion in the prior. Default is to set to TRUE}
\item{print_Z}{An option to output a summary of the allocations after each iteration. Default is to not print this summary.}
\item{num.cores}{Number of cores to be used by the function. Default is to use 4 cores.}
\item{mu_star_1_J_initial}{Initial estimates of mean expressions by bayNorm.}
\item{phi_star_1_J_initial}{Initial estimates of dispersions by bayNorm.}
}
\description{
Function to be used if want to obtain more posterior samples after running normHDP_mcmc.
}
\examples{
# Load data
load('data/case1_mcmc.RData')

# Run further MCMC
normHDP_mcmc_2(normHDP_output = case1_mcmc,
               number_iter = 5000)
}
