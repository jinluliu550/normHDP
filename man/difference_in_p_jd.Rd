\name{difference_in_p_jd}
\alias{difference_in_p_jd}
\title{Measure difference between component probabilities of each mouse pair}
\usage{
difference_in_p_jd(p_jd_trace,
                   alpha = 0.95)
}
\arguments{
\item{p_jd_trace}{trace of mouse-specfic component probabilities.}
\item{alpha}{threshold for classification, default is set to 0.95.}
}
\value{
\item{all_obs}{data frame of details for all mouse pairings and clusters.}
\item{significant_obs}{data frame of details for mouse pairings and clusters which are classified to show sigificant difference between mouse pairs.}
\item{plot}{a scatter plot to show the relationship between mean of absolute difference and tail probabilities.}
}
\description{
Measure difference between component probabilities of each mouse pairs using all posterior samples.
}
\examples{

# An example with threshold set to 0.95.

difference_in_p_jd(p_jd_trace = normHDP_post_MCMC$P_J_D_output,
                   alpha = 0.95)
}
