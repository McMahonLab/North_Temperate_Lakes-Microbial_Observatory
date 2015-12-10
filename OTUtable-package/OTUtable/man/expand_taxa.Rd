\name{expand_taxa}
\alias{expand_taxa}

\title{
Expand taxonomy assignment into multiple columns
}
\description{
The taxonomy file from mothur includes the assignment as a single string. This command takes that string and splits it into separate columns for each phylogenetic level. This makes searching for particular groups and combining OTUs based on taxonomy easier.
}
\usage{
expand_taxa(taxonomy_file)
}
\arguments{
	\item{taxonomy_file}{
  Path to a taxonomy file that has been processed using clean_taxonomy() and saved as a .csv file
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
#taxonomy <- expand_taxa(path)

}
