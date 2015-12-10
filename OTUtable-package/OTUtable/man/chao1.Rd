\name{chao1}
\alias{chao1}

\title{
Chao1 Richness
}
\description{
Calculates Chao1 richness of a vector of relative abundance data. This alpha diversity metric takes into account the number of singletons and doubletons for a more accurate estimate than observed richness.
}
\usage{
chao1(sample)
}
\arguments{
	\item{sample}{
  A vector of relative abundance data, typically a column in a matrix
}
}

\value{
Returns a single number indicating the estimated richness in the tested sample based on the number of taxa appearing only once or twice
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
Use apply functions to calculate Chao1 richness for all samples in a matrix

}

\keyword{ biodiversity }

\examples{
data(otu_table)

chao1_richness <- apply(otu_table, 2, chao1)
}
