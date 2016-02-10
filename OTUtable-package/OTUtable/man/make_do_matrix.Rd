\name{make_do_matrix}
\alias{make_do_matrix}

\title{
Make matrix of dissolved oxygen data
}
\description{
Takes a given sample ID and converts the dissolved oxygen data in data(metadata) from long format into a matrix. This is useful for plotting using plot_column()
}
\usage{
make_do_matrix(sampleID, field_data)
}
\arguments{
	\item{sampleID}{
  A regular expression used to select a group of samples
}
	\item{field_data}{
  A dataset of DO profiles in long format. Column names must be the same as the metadata file provided with this package
}
}
\details{
Also fills in NA values with the average of the depth above and below the missing value. If the value is at the bottom of the water column, the second deepest is substituted.
}
\value{
Returns matrix of DO data with depth in rows and date in columns
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
This is mainly used for generating contour plots. In general, long format is easier to work with.
In the metadata file included in this package, each DO measurement is listed twice, once under the epilimnion sample name and again under the hypolimnion sample name.
}

\keyword{water_column_plots}

\examples{
data(metadata)

dissolved_oxygen <- make_do_matrix("TBE.....07", metadata)
}
