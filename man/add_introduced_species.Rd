% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_introduced_species.R
\name{add_introduced_species}
\alias{add_introduced_species}
\title{add_introduced_species}
\usage{
add_introduced_species(
  native_parameters,
  introduced_lower_bound_growth_rate = 1,
  introduced_upper_bound_growth_rate = 5,
  introduced_self_sign,
  introduced_row_signs,
  introduced_col_signs,
  introduced_k
)
}
\arguments{
\item{native_parameters}{An object of type ´EEM´.}

\item{introduced_lower_bound_growth_rate}{Lower bound for the growth rate of
the introduced species.}

\item{introduced_upper_bound_growth_rate}{Upper bound for the growth rate of
the introduced species.}

\item{introduced_self_sign}{Desired sign of the self-interaction of the
introduced species.}

\item{introduced_row_signs}{A vector of signs for the interactions of the
introduced species with the native species (row interactions).}

\item{introduced_col_signs}{A vector of signs for the interactions of the
introduced species with the native species (column interactions).}

\item{introduced_k}{The equilibrium value for the introduced species.}
}
\value{
A list of extended parameters with the introduced species added.
}
\description{
This function adds an introduced species to a list of native species
parameters.
It samples the growth rate of the introduced species from a uniform
distribution
between specified lower and upper bounds, and samples the interactions
between
the introduced species and the native species, ensuring that the
self-interaction sign matches the specified \code{introduced_self_sign} and
a carrying capacity \code{k} is reached.
}
\details{
The function samples the growth rate of the introduced species from
a uniform distribution between the specified lower and upper bounds. It also
samples the interactions between the introduced species and the native
species, ensuring that the self-interaction sign matches the specified
\code{introduced_self_sign}. The interactions are added to the existing
\code{interaction_matrix} of each native species, and the growth rates are
extended to include the introduced species.
}
\examples{
library(EEMtoolbox)
output <- EEMtoolbox::EEM(matrix(c(-1, -1, 1, -1), ncol = 2), n_ensemble = 2)
intro <- add_introduced_species(output,
                                introduced_self_sign = -1,
                                introduced_row_signs = c(1, 1),
                                introduced_col_signs = c(-1, -1),
                                introduced_k = 10)
}
\author{
Maude Vernet
}
