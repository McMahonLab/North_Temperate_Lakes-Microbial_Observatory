\name{make_temp_matrix}
\alias{make_temp_matrix}

\title{
Make matrix of temperature data
}
\description{
Takes a given sample ID and converts temperature data of water profiles over time from long format into a matrix. This is most often useful for plotting using plot_column().
}
\usage{
make_temp_matrix(sampleID, field_data)
}
\arguments{
	\item{sampleID}{
  A regular expression used to select a group of samples
}
	\item{field_data}{
  A dataset of temperature profiles in long format. Column names must be the same as the metadata file provided with this package
}
}

\value{
Returns matrix of temperature data with depth in rows and date in columns
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
This is mainly used for generating contour plots. In general, long format is easier to work with.
In the included metadata file, each temperature measurement is recorded twice, once as epilimnion and once as hypolimnion.
}

\keyword{water_column_plots}

\examples{
data(metadata)

temp <- make_temp_matrix("TBE.....07", metadata)
}
