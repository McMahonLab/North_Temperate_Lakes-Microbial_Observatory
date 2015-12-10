\name{zscore}
\alias{zscore}

\title{
Z-score normalize relative abundance data
}
\description{
Normalizes taxa abundances in a table of relative abundance data using the z-score method.
}
\usage{
zscore(table)
}
\arguments{
	\item{table}{
  A table of relative abundance data with taxa in rows and samples in columns
}
}

\value{
Returns a table with relative abundance data replaced by z-scores
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
There is debate on whether this method of normalization is valid for microbial communities, as their abundance distrubtions tend to be heavily skewed.
The author found it useful for plotting heatmaps and for input into network analysis.
}

\keyword{data_processing}

\examples{
data(otu_table)
z_otu_table <- zscore(otu_table)
}
