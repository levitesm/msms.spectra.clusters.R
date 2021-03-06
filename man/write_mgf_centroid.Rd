% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{write_mgf_centroid}
\alias{write_mgf_centroid}
\title{A Write MGF Centroid function}
\usage{
write_mgf_centroid(mgf, path)
}
\arguments{
\item{mgf}{The mgf data.table object}

\item{path}{A path to the .mgf file}
}
\description{
This function allows you to write an .mgf data.table into an .mgf file. MGF should be of a Centroid type.
}
\examples{
write_mgf_centroid(mgf1,'C:/docs/mgf1.mgf')
}
\keyword{write}
