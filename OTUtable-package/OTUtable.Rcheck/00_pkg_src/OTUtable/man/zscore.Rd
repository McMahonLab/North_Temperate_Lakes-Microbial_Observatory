\name{zscore}
\alias{zscore}

\title{
Z-score normalize relative abundance data
}
\description{
Normalizes taxa abundances in a table of relative abundance data using the z-score method.
((Abundance of one OTU in one sample) - (mean abundance for that OTU ))/(standard deviation of that OTU)
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
I found it useful for plotting heatmaps and for input into network analysis.
}

\keyword{data_processing}

\examples{
data(otu_table)

# Create a small table for z-score normalization
example <- year_subset("05", otu_table)
example <- bog_subset("TBE", example)

# Remove OTUs that are not present in this subset
example <- example[which(rowSums(example) > 0), ]

z_otu_table <- zscore(example)
}
