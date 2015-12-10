\name{pielou}
\alias{pielou}

\title{
Pielou's Evenness
}
\description{
Calculates Pielou's evenness for a single vector of relative abundance data
}
\usage{
pielou(sample)
}

\arguments{
  \item{sample}{
  A vector of relative abundance data
}
}

\value{
  Returns a single value indicating the evenness of a community
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\keyword{biodiversity}
\examples{
data(otu_table)
even <- apply(otu_table, 2, pielou)

}


