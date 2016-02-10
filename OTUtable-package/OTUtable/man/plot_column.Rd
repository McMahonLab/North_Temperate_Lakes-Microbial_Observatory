\name{plot_column}
\alias{plot_column}

\title{
Plot DO or temperature data from a depth profile over time
}
\description{
Takes output from make_do_matrix or make_temp_matrix and plots using filled.contour()
}
\usage{
plot_column(data_matrix, title)
}
\arguments{
	\item{data_matrix}{
  A matrix output by make_do_matrix() or make_temp_matrix()
}
	\item{title}{
  The title you would like on the plot
}
}

\value{
Plots a filled contour plot showing the water column over time
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}
\note{
Depends on the function rotate(). The functions make_do_matrix() and make_temp_matrix() fill in missing values with the average of the measurement at each depth above and below; however, if missing values are present in the matrix for plotting, these will appear as white space on the plot.
}

\keyword{water_column_plots}

\examples{
data(metadata)
temp <- make_temp_matrix("TBE.....07", metadata)
plot_column(temp, "Trout Bog 2007 Temperature")
}
