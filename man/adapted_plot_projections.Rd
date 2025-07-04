% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapted_plot_projections.R
\name{adapted_plot_projections}
\alias{adapted_plot_projections}
\title{adapted_plot_projections}
\usage{
adapted_plot_projections(
  projections,
  title = NA,
  multiplier = 1,
  scaled = TRUE
)
}
\arguments{
\item{projections}{A data frame containing the projections of species
abundance.}

\item{title}{Title of the plot.}

\item{multiplier}{A numeric value to scale the abundance values, should be
the same as the one used when running \code{EEM}.}

\item{scaled}{Logical value indicating whether to scale the y-axis.}
}
\value{
A ggplot object with the projections of species abundance.
}
\description{
Plot projections of species abundance
}
\author{
Maude Vernet
}
