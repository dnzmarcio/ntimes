% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{helper_mann_whitney}
\alias{helper_mann_whitney}
\title{Helper for Mann-Whitney test}
\usage{
helper_mann_whitney(x, g, paired, alternative, conf.level)
}
\arguments{
\item{x}{a numeric vector.}

\item{g}{a factor.}

\item{paired}{a logical value.}

\item{alternative}{a character value.}

\item{conf.level}{a numeric value.}
}
\value{
a list with \code{test} indicating the test that was
performed, its respective \code{p.value} and \code{lower} and \code{upper}
limits of the confidence interval for the effect size considered in the test.
}
\description{
It performs the Mann-Whitney test.
}
\details{
This function can be modified by the user,
but input and output should be kept the same.
}
