\name{global_marker_genes}
\alias{global_marker_genes}
\title{Classifcation of global marker genes}
\usage{
global_marker_genes(normHDP_post_output,
                    threshold,
                    alpha_M = NULL,
                    alpha_D = NULL,
                    num.cores,
                    run.on.pc = TRUE)
}
\arguments{
\item{normHDP_post_output}{direct output from function normHDP_MCMC_fixed_z().}
\item{threshold}{a list of threshold for the computation of tail probabilities. The list has length 2, the first and the second item of the list is the threshold corresponding to the mean expression and dispersions, respectively.}
\item{alpha_M}{threshold for the classification of global differentially expressed genes. By default, the thereshold is set to control the expected false discovery rate (EFDR) to 5 percent.}
\item{alpha_D}{threshold for the classification of global differentially dispersed genes. By default, the thereshold is set to control the expected false discovery rate (EFDR) to 5 percent.}
\item{num.cores}{number of cores.}
\item{run.on.pc}{an option to run the algorithm on PC. Default is set to TRUE.}
}
\value{
\item{marker_DE}{summary statisics including the absolute log-fold change, tail probabilities and classification of global DE for each gene.}
\item{marker_DD}{summary statisics including the absolute log-fold change, tail probabilities and classification of global DD for each gene.}
\item{alpha_M}{threshold for the classification of global differentially expressed genes. By default, the thereshold is set to control the expected false discovery rate (EFDR) to 5 percent.}
\item{alpha_D}{threshold for the classification of global differentially dispersed genes. By default, the thereshold is set to control the expected false discovery rate (EFDR) to 5 percent.}
}
\description{
Finding the global marker genes
}
\examples{
global_marker_genes(normHDP_post_output = normHDP_post_MCMC,
                    threshold = c('mu' = 1.5, 'phi' = 1.5),
                    num.cores = 10)
}
