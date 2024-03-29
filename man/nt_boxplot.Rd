% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot.R
\name{nt_boxplot}
\alias{nt_boxplot}
\title{Boxplot}
\usage{
nt_boxplot(
  data,
  group = NULL,
  labels = NULL,
  save = FALSE,
  fig.height = 5,
  fig.width = 5,
  std_fun = std_boxplot,
  std_fun_group = std_boxplot_group,
  ...
)
}
\arguments{
\item{data}{a data frame with the variables.}

\item{group}{a character value indicating the group variable.}

\item{labels}{a list of labels with components given by their variable names.}

\item{save}{a logical value indicating whether the output
should be saved as a jpeg file.}

\item{fig.height}{a numeric value indicating the height (in) of the file.}

\item{fig.width}{a numeric value indicating the width (in) of the file.}

\item{std_fun}{a function to plot a boxplot when \code{group = NULL}.
It must follow the same structure of \code{\link{std_boxplot}}.}

\item{std_fun_group}{a function to plot a boxplot when \code{group}
is provided. It must follow the same structure of
\code{\link{std_boxplot_group}}.}

\item{...}{additional input arguments that may be used when creating your own function.}
}
\value{
a list of ggplot objects with each item named by the column names from
\code{var}.
}
\description{
Plot boxplot for several variables.
}
\details{
The functions \code{\link{std_boxplot}} and
\code{\link{std_boxplot_group}} can be modified by the user in order to
customize the boxplots a prior.
The plots also can be modified a posterior as a regular ggplot object.
See \code{\link[ggplot2]{geom_boxplot}}, \code{\link{std_boxplot}} and
\code{\link{std_boxplot_group}}.
}
\examples{
data(iris)

iris |> nt_boxplot(group = Species)

}
