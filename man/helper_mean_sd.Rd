% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{helper_mean_sd}
\alias{helper_mean_sd}
\title{Mean and standard deviation}
\usage{
helper_mean_sd(var, digits, ...)
}
\arguments{
\item{var}{a numeric vector.}

\item{digits}{a numeric value specifying the number of digits to present the results.}

\item{...}{additional input arguments that may be used when creating your own function.}
}
\value{
a list with the first element \code{name} as the measure name and the
second element as the \code{value} for a given variable.
}
\description{
It calculates mean and standard deviation, concatenating them to
present on a table.
}
\details{
This function can be modified by the user,
but input and output should be kept the same.
}
