\name{strat_metric}
\alias{strat_metric}

\title{
Plot temperature difference between 0.5 meters and max depth
}
\description{
Used as a metric of stratification. Calculates the difference in temperature just below the surface (0.5m) and the maximum depth.
}
\usage{
strat_metric(data_matrix)
}
\arguments{
	\item{data_matrix}{
  A matrix output by make_do_matrix or make_temp_matrix
}
	
}

\value{
Returns a vector of differences between the two depths
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{water_column_plots}

\examples{
data(metadata)
temp <- make_temp_matrix("TBE.....07", metadata)

diff <- strat_metric(temp)
}
