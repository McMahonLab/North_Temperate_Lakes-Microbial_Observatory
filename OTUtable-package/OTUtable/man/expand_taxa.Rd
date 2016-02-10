\name{expand_taxa}
\alias{expand_taxa}

\title{
Expand taxonomy assignment into multiple columns
}
\description{
The taxonomy file from mothur includes the assignment as a single string. This command takes that string and splits it into separate columns for each phylogenetic level. This makes searching for particular groups and combining OTUs based on taxonomy easier. There is also an option to remove the bootstrap values included in the full taxonomic assignments. These were removed in the NTL-Microbial Observatory taxonomy file included in this R package, but are still available in the data files included in the NTL-Microbial Observatory GitHub repository. 

The vocabulary for classifcation in the freshwater database is used for this command. In the output of expand_taxa(), "Lineage" = family level or 95 percent identity in the 16S region, "Clade" = genus level or 97 percent identity, and "Tribe" = species level, or 98 percent identity.
}
\usage{
expand_taxa(taxonomy_file, remove_bootstrap)
}
\arguments{
	\item{taxonomy_file}{
  Path to a taxonomy file that has been processed using clean_taxonomy() and saved as a .csv file
}
\item{remove_bootstrap}{
  TRUE or FALSE. If true, removes bootstrap values from the taxonomic assignments. This is useful when matching taxonomic assignments that would be identical except for their bootstrap values.
}
}

\value{
Returns a dataframe with each phylogenetic level in a separate column; each row indicates an OTU.
}


\author{
Alexandra Linz <amlinz16@gmail.com>
}

\keyword{format_mothur_output}

\examples{

#Example path only: path <- "mothur_output/bogs_taxonomy.csv"
#taxonomy <- expand_taxa(path, remove_bootstrap = T)

}
