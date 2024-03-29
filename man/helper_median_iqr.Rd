% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{helper_median_iqr}
\alias{helper_median_iqr}
\title{Median with first and third quantiles}
\usage{
helper_median_iqr(var, digits, ...)
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
It calculates median with quantiles 25% and 75%, concatenating them
to present on a table.
}
\details{
This function can be modified by the user,
but input and output should be kept the same.
}
