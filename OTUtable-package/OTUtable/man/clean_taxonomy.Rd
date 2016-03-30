\name{clean_taxonomy}
\alias{clean_taxonomy}

\title{
Clean Taxonomy File
}
\description{
Reduces information in a mothur .taxonomy file by removing the second column with the number of reads per OTU. It will also check for and remove OTUs in the taxonomy file that are not in OTU table - this may be the case if rarefaction was performed after classification, as was used in the NTL-Microbial Observatory dataset.
}
\usage{
clean_taxonomy(taxonomy_file, table, remove_bootstrap)
}
\arguments{
	\item{taxonomy_file}{
  A .taxonomy file output by the McMahon Lab 16S classification workflow
}
	\item{table}{
  An OTU table containing OTU numbers as row names
}
  \item{remove_bootstrap}{
  TRUE or FALSE: if TRUE, removes bootstrap values from the classification strings
  }
}

\value{
Returns the taxonomy with OTUs as row names and seven columns containing each taxonomic level (Kingdom, Phylum, Class, Order, Lineage, Clade, and Tribe)
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{ format_mothur_output }

\examples{

# Example path only: path <- "mothur_output/bogs.taxonomy"
# table <- clean_shared("mothur_output/bogs.shared", trim.names = T)
# taxonomy <- clean_taxonomy(path, table, remove_bootstrap = F)

}

