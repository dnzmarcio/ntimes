% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kaplan_meier.R
\name{nt_km}
\alias{nt_km}
\title{Kaplan-Meier plot}
\usage{
nt_km(
  data,
  time,
  status,
  labels = NULL,
  xlab = "Time",
  ylab = "Survival",
  risktable.title = "n at risk",
  save = FALSE,
  fig.height = 5,
  fig.width = 5,
  std_fun = std_km,
  std_fun_group = std_km_group,
  time.points = NULL,
  format = TRUE,
  digits = 2,
  file = "survival"
)
}
\arguments{
\item{data}{a data frame with the variables.}

\item{time}{a data frame with the time to event variable.}

\item{status}{a data frame with the indicator associated to events.}

\item{labels}{a list of labels with components given by their variable names.}

\item{xlab}{a character value specifying the x axis label.}

\item{ylab}{a character value specifying the y axis label.}

\item{save}{a logical value indicating whether the output
should be saved as a jpeg file.}

\item{fig.height}{a numeric value indicating the height (in) of the file.}

\item{fig.width}{a numeric value indicating the width (in) of the file.}

\item{std_fun}{a function to plot a barplot when \code{group = NULL}.
It must follow the same structure of \code{\link{std_barplot}}.}

\item{std_fun_group}{a function to plot a dotplot when \code{group}
is provided. It must follow the same structure of
\code{\link{std_barplot_group}}.}

\item{time.points}{a numeric vector of time points to evaluate the survival curves.}

\item{format}{a logical value indicating whether the output should be formatted.}

\item{digits}{a numerical value defining of digits to present the results.}

\item{file}{a character indicating the name of output file in csv format to be saved.}
}
\value{
a list of ggplot objects with each item named by the column names from
\code{var}.
}
\description{
Plot Kaplan-Meier curves for several variables.
}
\details{
The functions \code{\link{std_km}} and
\code{\link{std_km_group}} are standard functions that can be
modified by the user in order to customize the barplots a prior.
The plots also can be modified a posterior as a regular ggplot object.
See \code{\link{std_km}} and \code{\link{std_km_group}}.
}
\examples{
library(survival)
data(lung)

lung_nt <- lung \%>\% mutate(sex = factor(sex, levels = 1:2,
                                    labels = c("Female", "Male")),
                       ph.ecog = as.factor(ph.ecog)) \%>\%
                       select(sex, ph.ecog, time, status)
lung_nt \%>\% nt_km(time = time, status = status,
                 labels = list(sex = "Sex", ph.ecog = "ECOG"))

}
