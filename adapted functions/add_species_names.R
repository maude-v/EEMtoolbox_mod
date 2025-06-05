#' @title add_species_names
#' @description
#' Add species names to a numeric vector or matrix
#' This helper takes an ´EEM´ object and adds species names to the
#' interaction matrix.
#' @param parameter a set of parameters obtained with {EEMtoolbox::EEM()}.
#' @param species_names Character vector of length p. The new names to assign.
#' @return The same object \code{parameter} but with the names set to
#' \code{species_names} on the interaction matrix.
#' @examples
#' dingo_EEM <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
#' add_species_names(dingo_matrix,
#'                   c("sp1", "sp2", "sp3", "sp4", "sp5", "sp6", sp7", "sp8"))
#' @author Maude Vernet
#' @export
add_species_names <- function(parameter,
                              species_names) {
  lapply(parameter, function(x) {
    colnames(x$interaction_matrix) <- c(species_names)
    rownames(x$interaction_matrix) <- c(species_names)
    return(x)
  })
}
