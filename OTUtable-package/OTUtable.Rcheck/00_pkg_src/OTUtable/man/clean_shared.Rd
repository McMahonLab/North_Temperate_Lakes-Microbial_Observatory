\name{clean_shared}
\alias{clean_shared}

\title{
Reformat a shared file
}
\description{
Converts a mothur .shared file into a simplified OTU table. The columns indicating total reads for each OTU and the clustering level are removed, and the table is transposed so that OTUs are rows and samples are columns. The "trim.names" variable provides an option to shorten sample names to the first "." character - this is specific to the NTL-Microbial Observatory dataset. Manual curation of sample names took place after this step for the NTL-Microbial Observatory dataset in order to maintain consistency across all sample names.
}
\usage{
clean_shared(shared_file, trim.names)
}
\arguments{
	\item{shared_file}{
  A .shared file output by mothur
}
\item{trim.names}{
  TRUE or FALSE - if TRUE, sample names will be trimmed to the first "." character.
}
}

\value{
Returns an OTU table with samples as columns and OTUs as rows. 
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{ format_mothur_output }

\examples{

# Example path only: path <- "mothur_output/bogs.shared"

# otu_table <- clean_shared(path, trim.names = T)
# write.csv(otu_table, file = "bogs_otu_table.csv", quote = F, row.names = T)

}
