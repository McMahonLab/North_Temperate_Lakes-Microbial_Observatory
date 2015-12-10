\name{extract_date}
\alias{extract_date}

\title{
Extract sampling date from a vector of sample names
}
\description{
The date each sample was collected is encoded in the sample ID. Extract this into R date format using this command. 
}
\usage{
extract_date(sample_ids)
}
\arguments{
	\item{sample_ids}{
  A vector of sample names. Samples must be labeled using the bog, layer, date, and replicate system (MAH04JUL05.R1 = Mary Lake Hypolimnion, 04Jul05, replicate 1)
}
}

\value{
Returns a vector of dates corresponding to each sample
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}


\keyword{data_processing}

\examples{

samples <- c("TBE01JUN09.R1", "TBE05JUN09", "TBE10JUN09.R2")
extract_date(samples)

}
