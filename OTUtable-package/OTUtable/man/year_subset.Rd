\name{year_subset}
\alias{year_subset}

\title{
Subset samples by a specific year
}
\description{
Takes the year value in the last two digits of the sample ID and allows selection of a single year of data.
}
\usage{
year_subset(year_id, table)
}
\arguments{
	\item{table}{
  A table containing the relative abundances of each taxa
}
\item{year_id}{
  Two digit code indicating the last two digits of the year of interest (05, 07, 08, 09) surrounded by quotes
}
}

\value{
Returns an OTU table containing only samples from the specified year
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{

Regular expressions are permitted.
}

\keyword{data_processing}

\examples{
data(otu_table)

seven <- year_subset("07", otu_table)
}
