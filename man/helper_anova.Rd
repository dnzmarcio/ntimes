% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{helper_anova}
\alias{helper_anova}
\title{Helper for Analysis of Variance}
\usage{
helper_anova(x, g)
}
\arguments{
\item{x}{a numeric vector.}

\item{g}{a factor.}
}
\value{
a list with \code{test} indicating the test that was
performed and its respective \code{p.value}.
}
\description{
It performs the Analysis of Variance.
}
\details{
This function can be modified by the user,
but input and output should be kept the same.
}
