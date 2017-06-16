\name{rotate}
\alias{rotate}

\title{
Rotate a matrix
}
\description{
Rotates a matrix of data so that columns are reversed
}
\usage{
rotate(data_matrix)
}
\arguments{
	\item{data_matrix}{
  Used in this package with matrix output by make_do_matrix or make_temp_matrix as part of the function plot_column(). Any matrix will work, though.
}
	
}
\details{
Used to rotate the DO or temperature matrices so that depth 0 is at the top of a contour plot and the max depth is at the bottom.
}
\value{
Returns a matrix that has been rotated so that it reads from bottom to top
}

\author{
An anonymous author on Stack Overflow
Alexandra Linz <amlinz16@gmail.com>
}
\note{
Used with make_do_matrix(), make_temp_matrix(), and plot_column(). plot_column() depends on this function.
}

\keyword{water_column_plots}

\examples{
data(metadata)
temp <- make_temp_matrix("TBE.....07", metadata)
r_temp <- rotate(temp)

}
