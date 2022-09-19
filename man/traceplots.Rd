\name{traceplots}
\alias{traceplots}
\title{Trace plot of parameters}
\usage{
traceplots(normHDP_mcmc_output,
           normHDP_post_mcmc_output,
           parameter,
           j = NULL,
           g = NULL,
           c = NULL,
           d = NULL)
}
\arguments{
\item{normHDP_mcmc_output}{Output from function normHDP_mcmc.}
\item{normHDP_post_mcmc_output}{Output from function normHDP_mcmc_post.}
\item{parameter}{parameter of interest. This should be one of 'mu' (mean expression - need to input j and g), 'phi' (dispersion - need to input j and g), 'beta' (capture effiencicies - need to input d and c), 'alpha', 'alpha_zero', 'b' and 'alpha_phi_2'.}
}
\value{
Trace plot of specified parameter.
}
\description{
To draw a trace plot of parameter of interest.
}
\examples{
# Load data
load('data/case1_mcmc.RData')
load('data/case1_mcmc_post.RData)

# Trace plot of mu
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'mu',
           j = 10,
           g = 10)

# Trace plot of phi
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'phi',
           j = 10,
           g = 10)

# Trace plot of alpha
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'alpha')

# Trace plot of alpha
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'alpha_zero')

# Trace plot of beta
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'beta',
           d = 1,
           c = 10)

# Trace plot of alpha_phi_2
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'alpha_phi_2')

# Trace plot of b
traceplots(normHDP_mcmc_output = case1_mcmc,
           normHDP_post_mcmc_output = case1_mcmc_post,
           parameter = 'b')
}
