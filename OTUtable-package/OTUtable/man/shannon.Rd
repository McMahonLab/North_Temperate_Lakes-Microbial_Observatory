\name{shannon}
\alias{shannon}

\title{
Shannon's Biodiversity Index
}
\description{
Calculates Shannon's Biodiversity Index on a single column of relative abundance data. This metric takes into account both richness and evenness.
}
\usage{
shannon(sample)
}
\arguments{
	\item{sample}{
  A vector of relative abundance data, typically a single column in a matrix
}
}

\value{
Returns a single number indicating the amount of biodiversity in the tested sample
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
Use apply functions to calculate Shannon's index for all samples in a matrix

}

\keyword{ biodiversity }

\examples{
data(otu_table)
richness <- apply(otu_table, 2, shannon)
}
