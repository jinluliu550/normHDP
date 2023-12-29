\name{ppp_calculation}
\alias{ppp_calculation}
\title{Computation of p-values for posterior predictive checks}
\usage{
ppp_calculation(normHDP_post_output,
                Y,
                num.cores = 1,
                rep.number,
                type = c('standard', 'mixed'),
                run.on.pc = TRUE)
}
\arguments{
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{Y}{input data.}
\item{num.cores}{number of cores to use.}
\item{rep.number}{total number of replicated datasets.}
\item{type}{method used to simulate the replicated data. If type is set to 'standard', then the dispersion parameters are drawn directly from the list of posterior samples to simulate the replicated data. If type is set to 'mixed', then the dispersion parmaters are simulated from the posterior samples of mean expressions and hyper-parameters following a log-normal distribution, and simulate the replicated data using the replicated dispersion parameters.}
}
\item{run.on.pc}{an option to choose where to run the algorithm. Important for parallel computing. By default, this is set to TRUE.}
\value{
Computation of p-values of posterior predictive checks. There are three discrepancy being measured: 1) the median difference between the observed and expected Y for each gene, 2) the IQR of the difference between the observed and expected Y for each gene, and 3) expected difference between the true and expected dropout probabilities for each gene.
}
\description{
The function computes p-values of posterior predictive checks with multiple replicates.
}
\examples{
# An example with 100 replicated datasets, running on 10 cores on a pc.

ppp_calculation(normHDP_post_output = normHDP_post_mcmc,
                Y = Y,
                num.cores = 10,
                rep.number = 100,
                type = 'mixed',
                run.on.pc = TRUE)
}
