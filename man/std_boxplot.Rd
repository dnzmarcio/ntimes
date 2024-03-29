% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot.R
\name{std_boxplot}
\alias{std_boxplot}
\title{Standard boxplot}
\usage{
std_boxplot(var, var.label, ...)
}
\arguments{
\item{var}{a numeric vector.}

\item{var.label}{a character value specifying the variable label.}

\item{...}{additional input arguments that may be used when creating your own function.}
}
\value{
a ggplot object.
}
\description{
A function to plot a boxplot without groups.
}
\details{
This function defines the standard boxplot without groups to be
plotted by the function \code{\link{nt_boxplot}}. It can be modified by the
user. See more details in \code{\link[ggplot2]{geom_boxplot}}.
}
