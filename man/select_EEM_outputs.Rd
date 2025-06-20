% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_EEM_outputs.R
\name{select_EEM_outputs}
\alias{select_EEM_outputs}
\title{select_EEM_outputs}
\usage{
select_EEM_outputs(
  ensemble,
  target_lower,
  target_upper,
  mode = "native",
  n_intro = 1
)
}
\arguments{
\item{ensemble}{an object of type \code{EEM}.}

\item{target_lower}{A numeric vector (or scalar) specifying the minimum
acceptable equilibrium abundance for each species.}

\item{target_upper}{A numeric vector (or scalar) specifying the maximum
acceptable equilibrium abundance for each species.}

\item{mode}{A character string indicating the mode of the simulation, either
\code{"native"} or \code{"disturbed"}. This affects how the equilibrium is evaluated.}

\item{n_intro}{An integer specifying the number of introduced species in the
simulation. This is used to label species in the output plot.}
}
\value{
A list containing two elements:
\itemize{
\item \code{outputs_selected}: Indices of the parameter sets that meet the criteria.
\item \code{outputs_unselected}: Indices of the parameter sets that do not meet the
criteria.
}
}
\description{
This function filters a list of ecological ensemble model (EEM) outputs
based on whether the system has reached equilibrium and if the
equilibrium abundances of species fall within specified bounds.
}
\details{
This function examines each element of \code{outputs} to (1) verify that the
system has reached equilibrium, and (2) check that every species’ equilibrium
abundance lies between \code{lower_bound} and \code{upper_bound}.
}
\examples{
output <- EEMtoolbox::EEM(matrix(c(-1, -1, 1, -1), ncol = 2), n_ensemble = 2)
output <- select_EEM_outputs(output,
                             target_lower  = rep(0.1, 2),
                             target_upper  = rep(15, 2))
}
\author{
Maude Vernet
}
