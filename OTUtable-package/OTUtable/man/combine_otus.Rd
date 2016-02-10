\name{combine_otus}
\alias{combine_otus}

\title{
Combine OTUs based on identical taxonomic assignments
}
\description{
Sums the abundances of OTUs with the same taxonomy at a given level into a single vector for that taxonomy. This creates a new table of relative abundance data at a higher taxonomic level than OTU. This function only works with the OTU level as input, but can be used on any subset of the OTU table created by year_subset() or bog_subset().
If bootstrap values were not removed by expand_taxa(), this command will likely create spurious groupings based on identical bootstrap values.
}
\usage{
combine_otus(level, table, taxonomy)
}
\arguments{
	\item{level}{
  The desired level at which to combine OTUs; options are the column names from the taxonomy dataset
}
	\item{table}{
  An OTU table containing the relative abundances of each OTU.
}
	\item{taxonomy}{
  A taxonomy dataset in the form produced by expand_taxa().
}
}

\value{
Returns a table of relative abundance data with each row representing all OTUs of a given taxonomic assignment summed together. Row names are now the full taxonomic assignment of each row. To keep only the the lowest taxonomic level in the row names, run the function reduce_names()
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{data_processing}

\examples{

data(otu_table)
data(taxonomy)

example_table <- year_subset("05", otu_table)
example_table <- bog_subset("TBE", example_table)

phylum_table <- combine_otus("Phylum", example_table, taxonomy)

}
