\name{reduce_names}
\alias{reduce_names}

\title{
Shorten taxonomic assignment in table row names
}
\description{
Reduces the full string indicating taxonomy to the last classified level. Works on tables at levels higher than OTUs.
}
\usage{
reduce_names(table)
}
\arguments{
	\item{table}{
  A table containing the relative abundances of each taxa produced by combine_otus()
}
}

\value{
Returns the same table with shortened row names
}
\note{
This function is often most useful for plotting, so that the full string does not appear on the plot
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{data_processing}

\examples{
data(otu_table)
data(taxonomy)

# Create a small table for the example
example <- year_subset("05", otu_table)
example <- bog_subset("TBE", example)

clade_table <- combine_otus("Clade", example, taxonomy)
clade_table <- clade_table[which(rowSums(clade_table) > 0),]
head(rownames(clade_table))
reduced_clades <- reduce_names(clade_table)
head(rownames(reduced_clades))
}
