\name{clean_shared}
\alias{clean_shared}

\title{
Reformat mothur shared file
}
\description{
Converts a mothur .shared file into a simplified OTU table
}
\usage{
clean_shared(shared_file)
}
\arguments{
	\item{shared_file}{
  A .shared file output by mothur
}
}

\value{
Returns the shared file as an OTU table with samples as columns and OTUs as rows. Removes information about number of reads per sample and clustering level.
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{ format_mothur_output }

\examples{

#Example path only: path <- "mothur_output/bogs.shared"

#otu_table <- clean_shared(path)
#write.csv(otu_table, file = "bogs_otu_table.csv", quote = F, row.names=T)

}
