\name{latent_counts_plot}
\alias{latent_counts_plot}
\title{Plot of estimated latent count}
\usage{
latent_counts_plot(latent_counts_output,
                   width.plot = 1000,
                   height.plot = 400,
                   file.name.1,
                   file.name.2 = NULL)
}
\arguments{
\item{latent_counts_output}{direct output from function latent_counts().}
\item{width.plot}{width of the heat-map.}
\item{height.plot}{height of the heat-map.}
\item{file.name.1}{name of the first heat-map.}
\item{file.name.2}{name of the second heat-map. In the case when the genes are not reordered, file.name.2 is set to NULL.}
}
\value{
In the case when the genes are reordered by the tail probabilities, the function outputs 2 plots. Plot 1 is the heat-map with cells reordered by the clustering estimate, and genes reordered by tail probabilities of mean expressions. Plot 2 is the heat-map with cells reordered by the clustering estimate, and genes reordered by tail probabilities of dispersions. In the case when the genes are not reordered, the function only outputs one heat-map of latent count with reordered cells.
}
\description{
The function plots heat-map of estimated latent counts.
}
\examples{
latent_counts_plot(latent_counts_output = latent_counts_df,
                   width.plot = 1000,
                   height.plot = 400,
                   file.name.1 = 'plot1.png',
                   file.name.2 = 'plot2.png')
}
