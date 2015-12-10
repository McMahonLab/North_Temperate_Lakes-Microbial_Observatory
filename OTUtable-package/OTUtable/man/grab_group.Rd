\name{grab_group}
\alias{grab_group}

\title{
Subset OTU table by taxonomy
}
\description{
Returns a table containing only taxa from a given phylogenetic group
}
\usage{
grab_group(group, level, table, taxonomy)
}
\arguments{
	\item{group}{
  The phylogenetic classification of interest (can be a regular expression)
}
	\item{level}{
  The phylogenetic level of the group of interest (must be a column name in the taxonomy file)
}
	\item{table}{
  A table containing the relative abundances of each OTU
}
	\item{taxonomy}{
  A taxonomy dataset in the form produced by expand_taxa()
}
}

\value{
Returns a table with all taxa of a given taxonomic assignment
}
\note{
This function must be run on the OTU level table. However, the output of this function can be run through combine_otus() to create a higher level table of results.
}
\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{data_processing}

\examples{
 
data(otu_table)
data(taxonomy)

acI <- grab_group("acI", "Clade", otu_table, taxonomy)
verruco <- grab_group("Verrucomicrobia", "Phylum", otu_table, taxonomy)
}
