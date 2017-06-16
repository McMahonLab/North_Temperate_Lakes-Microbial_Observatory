\name{obs_richness}
\alias{obs_richness}

\title{
Observed Richness
}
\description{
Calculates observed richness on a single column of relative abundance data. 
}
\usage{
obs_richness(sample)
}
\arguments{
	\item{sample}{
  A vector of relative abundance data, typically a single column in a matrix
}
}

\value{
Returns a single number indicating the number of taxa in the tested sample
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
Use apply functions to calculate richness for all samples in a matrix

}

\keyword{ biodiversity }

\examples{
 data(otu_table)
 richness <- apply(otu_table, 2, obs_richness)
}
