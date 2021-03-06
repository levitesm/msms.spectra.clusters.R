% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{Spec_list_corr}
\alias{Spec_list_corr}
\title{A Parse MGF function}
\usage{
Spec_list_corr(mass_list1, mass_list2, D = 0.04, mode = 1)
}
\arguments{
\item{mass_list1}{The first list of masses.}

\item{mass_list2}{The secound list of masses.}

\item{D}{A maximal difference between two masses for them to be considered the same. Delaults to 0.04}

\item{mode}{If mode=1 the correlation score will be returned, else the number of correlationg peaks is returned. The correlation score is the percent of correlating peaks out of the mean length of two lists. Defaluts to 1.}
}
\description{
This function calculates the correlation between two lists of peak masses.
}
\examples{
corr <- Spec_list_corr(list1, list2)
}
\keyword{correlation}
