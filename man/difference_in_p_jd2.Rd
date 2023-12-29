\name{difference_in_p_jd2}
\alias{difference_in_p_jd2}
\title{Measure difference between component probabilities of two datasets}
\usage{
difference_in_p_jd2(p_jd_trace,
                    alpha = 0.95)
}
\arguments{
\item{p_jd_trace}{trace of mouse-specfic component probabilities.}
\item{alpha}{threshold for classification, default is set to 0.95.}
}
\value{
\item{all_obs}{data frame of details for all clusters.}
\item{significant_obs}{data frame of details for clusters which are classified to show sigificant difference between mouse pairs.}
\item{plot}{a scatter plot to show the relationship between mean of absolute difference and probability of having the component probability of data 1 greater than data 2.}
}
\description{
Measure difference between component probabilities of two data using all posterior samples. Note the function could only be used in cases with 2 datasets. The clusters are classified into 3 types: over-represented for data 2, under-represented for data 2 and stable.
}
\examples{

# An example with threshold set to 0.95.

difference_in_p_jd2(p_jd_trace = normHDP_post_MCMC$P_J_D_output,
                    alpha = 0.95)
}
