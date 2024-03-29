\name{normHDP_mcmc}
\alias{normHDP_mcmc}
\title{MCMC Simulation of unknown parameters}
\description{
Generate posterior samples of parameters of interests.
}
\usage{
normHDP_mcmc(Y,
             J,
             number_iter,
             thinning = 5,
             empirical = TRUE,
             burn_in = 3000,
             quadratic = FALSE,
             iter_update = 100,
             beta.mean = 0.06,
             alpha_mu_2 = NULL,
             adaptive_prop = 0.1,
             BB_SIZE = TRUE,
             print_Z = FALSE,
             num.cores = 4,
             save.only.z = FALSE,
             baynorm.beta = NULL,
             run.on.pc = TRUE)
}
\arguments{
\item{Y}{input dataset of length d.}
\item{J}{total number of components.}
\item{number_iter}{total number of iterations to run.}
\item{thinning}{an option to save a subset of the samples.}
\item{burn_in}{the first number of iterations to ignore.}
\item{quadratic}{prior relationship between the unique parameters on the log-scale, default is to set to TRUE.}
\item{iter_update}{print text after a given number of iterations to monitor progress.}
\item{beta.mean}{prior mean of capture effiencies, default is set to 0.06.}
\item{alpha_mu_2}{prior variance of mean expressions.}
\item{adaptive_prop}{added extra variation in adaptive MH, default is set to 0.1.}
\item{BB_SIZE}{method used in the baynorm() function to estimate dispersions, default is set to TRUE.}
\item{print_Z}{an option to print out a summary of allocations after each iteration.}
\item{num.cores}{number of processing cores.}
\item{save.only.z}{an option to only save the estimated allocations.}
\item{baynorm.beta}{an option to manually input the capture effiencices for empirical estimations.}
\item{run.on.pc}{an option to run the algorthm on PC. Default is set to TRUE.}
}
\value{
The function returns a list of the following items:
\item{b_output}{samples of mean of regression slopes.}
\item{alpha_phi_2_output}{samples of variance of dispersions.}
\item{Z_output}{samples of allocations.}
\item{P_J_D_output}{samples of dataset specific component probabilities.}
\item{P_output}{samples of component probabilities.}
\item{alpha_output}{samples of concentration parameters to control difference in component probabilities between datasets.}
\item{alpha_zero_output}{samples of concentration parameterst to control difference in component probabilities between components.}
\item{mu_star_1_J_output}{samples of estimated mean expressions.}
\item{phi_star_1_J_output}{sampes of estimated dispersions.}
\item{Beta_output}{samples of capture efficiencies.}
\item{acceptance_prob_list}{a data frame to summarize acceptance probabilities of MH.}
\item{G}{number of genes.}
\item{D}{number of datasets.}
\item{C}{number of cells in each dataset.}
}
\examples{

# An example run with 3 components.

normHDP_mcmc(Y = Y,
             J = 3,
             number_iter = 3000,
             thinning = 5,
             empirical = TRUE,
             burn_in = 1000,
             quadratic = FALSE,
             iter_update = 100,
             beta.mean = 0.72,
             alpha_mu_2 = 20,
             print_Z = TRUE,
             num.cores = 4)

}
