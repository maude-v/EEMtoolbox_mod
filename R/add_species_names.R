#' @title add_species_names
#' @author Maude Vernet
#' @description
#' Add species names to a numeric vector or matrix
#' This helper takes an ´EEM´ object and adds species names to the
#' interaction matrix.
#' @param parameter a set of parameters obtained with \code{EEMtoolbox::EEM()}.
#' @param species_names Character vector of length p. The new names to assign.
#' @return The same object \code{parameter} but with the names set to
#' \code{species_names} on the interaction matrix.
#' @examples
#' output <- EEMtoolbox::EEM(matrix(c(-1, -1, 1, -1), ncol = 2), n_ensemble = 2)
#' withnames <- add_species_names(output,
#'                                c("sp1", "sp2"))
#' @export
add_species_names <- function(parameter,
                              species_names) {
  lapply(parameter, function(x) {
    colnames(x$interaction_matrix) <- c(species_names)
    rownames(x$interaction_matrix) <- c(species_names)
    return(x)
  })
}
