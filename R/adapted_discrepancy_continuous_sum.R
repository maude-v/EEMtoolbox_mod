#' @title Adapted disc_func for EEM
#' @author Maude Vernet
#' @description A wrapper around \code{EEMtoolbox::dicrepancy_continuous_sum()} that
#’ allows you to enforce upper and lower bounds on equilibrium species abundances.
#' @inheritParams discrepancy_continuous_sum
#' @param target_lower Numeric vector of length p. If provided, these values will
#'   serve as the minimum allowable equilibrium abundances for each of the p species.
#' @param target_upper Numeric vector of length p. If provided, these values will
#'   serve as the maximum allowable equilibrium abundances for each of the p species.
#'
#' @details
#' This is exactly the same logic as \code{dicrepancy_continuous_sum()} in EEMtoolbox,
#' except that the standard EEMtoolbox version does _not_ check for user‐specified
#' equilibrium bounds.
#'
#' `important notice`: numbers in `target_lower` and `target_upper` are to be
#' inputed directly, and not inside a variable. Otherwise, they are not found by
#' the function.
#'
#' @examples
#' out <- EEMtoolbox::EEM(
#'              interaction_matrix = matrix(c(-1, -1, 1, -1), ncol = 2),
#'              n_ensemble = 2,
#'              disc_func = function(data) {
#'              EEMtoolbox::adapted_discrepancy_continuous_sum(
#'              data,
#'              target_lower = rep(1,2),
#'              target_upper = rep(20,2))})
#'
#' @return A list with the same elements as \code{EEMtoolbox::dicrepancy_continuous_sum()}
#'
#' @seealso
#' \code{EEMtoolbox::dicrepancy_continuous_sum()} (original function)
#'
#' @export
adapted_discrepancy_continuous_sum <- function(data,
                                               target_lower = NULL,
                                               target_upper = NULL) {
  n_species <- length(data) / 2
  equilibrium_points <- data[seq_len(n_species)]
  # computed equilibrium abundances
  stability_eigenvalues <- data[(n_species + 1):length(data)]
  # stability eigenvalues

  # Feasibility penalty: sum of negative equilibrium parts.
  measured_infeasibility <- abs(sum(pmin(0, equilibrium_points)))
  # Stability penalty: sum of positive eigenvalues.
  measured_instability <- sum(pmax(0, stability_eigenvalues))
  #until here, this is the same function as the original one


  # Initialize target penalty. is 0 if no target equilibrium is provided
  target_penalty <- 0

  # If target lower and upper bounds are provided, calculate per-species penalties.
  if (!is.null(target_lower) && !is.null(target_upper)) {
    if (length(target_lower) != n_species || length(target_upper) != n_species)
    {
      stop("target_lower and target_upper must be vectors of
           length equal to the number of species.")
    }

    # For each species, if the computed equilibrium is below target_lower, add (target_lower - computed).
    lower_penalty <- pmax(0, target_lower - equilibrium_points)
    # For each species, if the computed equilibrium is above target_upper, add (computed - target_upper).
    upper_penalty <- pmax(0, equilibrium_points - target_upper)

    # Sum the penalties across all species.
    target_penalty <- sum(lower_penalty + upper_penalty)
  }

  # Total discrepancy is the sum of feasibility, instability, and target penalties.
  summ <- measured_infeasibility + measured_instability + target_penalty
  return(summ)
}
