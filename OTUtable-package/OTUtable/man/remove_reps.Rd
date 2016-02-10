\name{remove_reps}
\alias{remove_reps}

\title{
Remove the second replicate of each sample, when it exists
}
\description{
Sometimes it is desirable to remove replicate samples (often for plotting). This command removes all samples marked as replicate 2. Please note that you should always check the similarity of replicates for your metric of interest before removing them for aesthetic purposes.
}
\usage{
remove_reps(table)
}
\arguments{
	\item{table}{
  An OTU table containing the relative abundances of each OTU
}
}

\value{
Returns an OTU table containing only one replicate for each sample
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{data_processing}

\examples{
 data(otu_table)
 no_reps <- remove_reps(otu_table)
}
