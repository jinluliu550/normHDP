\name{cluster_summary}
\alias{cluster_summary}
\title{Plots to summarise point estimates of allocations}
\usage{
cluster_summary(Z_point)
}
\arguments{
\item{Z_point}{Point estimate of allocations.}
}
\value{
The output include 2 plots: Plot 1 is a bar chart summarizing the total number of cells allocated to each identified cluster. Plot 2 shows the proportion of cells from the first and second dataset for each cluster, and compare this proportion to the overall true proportion of cells.
}
\description{
To visualize point estimates of allocations.
}
\examples{
# Load data
load('data/case1_Z_estimate.RData')

# Obtain a summary of posterior allocations
cluster_summary(Z_point = case1_Z_estimate)
}
