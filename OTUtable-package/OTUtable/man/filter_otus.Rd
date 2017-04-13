\name{filter_taxa}
\alias{filter_taxa}

\title{
Filter Taxa Based on Abundance and Persistence
}
\description{
Returns a table containing only taxa that meet the imposed requirements of a minimum abundance and a minimum number of samples containing that taxon
}
\usage{
filter_taxa(table, abundance, persistence)
}
\arguments{
	\item{table}{
  A table containing the relative abundances of each OTU or taxon in the form produced by clean_shared(). Can be used on the output of grab_groups() or combine_otus()
}
	\item{abundance}{
  The minimum threshold for percentage of reads attributed to a taxon in at least one sample. Taxa at abundances greater than or equal this number will be retained.
}
	\item{persistence}{
  The minimum threshold for the percentage of samples in which a taxon has been observed.  Taxa at abundances greater than or equal this number will be retained.
}
}

\value{
Returns a table with all taxa that met the imposed thresholds
}
\note{
Thanks Juliana Dias for suggesting this function!

}
\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{data_processing}

\examples{

# To make a table containing only OTUs with at least 0.1% abundance 
# in at least one sample that were observed 
# (at any abundance) in at least 50% of samples:
# library(OTUtable)
# data(otu_table)
# filtered_table <- filter_taxa(otu_table, abundance = 0.1, persistence = 50)

# To make a table containing only phyla with at least 10% abundance 
# in any one sample and were observed 
# at any abundance in at least 10% of samples:
# data(taxonomy)
# phylum_table <- combine_otus("Phylum", otu_table, taxonomy)
# filtered_phylum_table <- filter_taxa(phylum_table, abundance = 10, persistence = 10)
}
