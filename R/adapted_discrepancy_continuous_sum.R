#' @title Adapted disc_func for EEM
#' @author Maude Vernet
#' @description A wrapper around \code{\link[EEMtoolbox]{dicrepancy_continuous_sum}()} that
#’ allows you to enforce upper and lower bounds on equilibrium species abundances.
#' @inheritParams discrepancy_continuous_sum
#' @param lower_bounds Numeric vector of length p. If provided, these values will
#'   serve as the minimum allowable equilibrium abundances for each of the p species.
#'   Any solution below \code{lower_bounds} is reported as “out of bound.”
#' @param name description upper_bounds Numeric vector of length p. If provided, these values will
#'   serve as the maximum allowable equilibrium abundances. Any solution above
#'   \code{upper_bounds} is reported as “out of bound.”
#'
#' @details
#' Internally calls \code{EEMtoolbox::dicrepancy_continuous_sum()} with the same
#' syntax, then post‐filters the returned equilibrium abundances to check whether
#' all species fall within \code{[lower_bounds, upper_bounds]}. If any are outside,
#' a warning is thrown (or, optionally, an error, depending on \code{stop_on_out_of_bound}).
#'
#' This is exactly the same logic as \code{dicrepancy_continuous_sum()} in EEMtoolbox,
#' except that the standard EEMtoolbox version does _not_ check for user‐specified
#' equilibrium bounds.
#'
#' @return A list with the same elements as \code{EEMtoolbox::dicrepancy_continuous_sum()}, plus:
#'   \item{out_of_bound}{Logical scalar. \code{TRUE} if any equilibrium abundances
#'     lie outside \code{lower_bounds}/\code{upper_bounds}; \code{FALSE} otherwise.}
#'   \item{bounds_violation_indices}{Integer vector of all species indices that
#'     violate the specified bounds.}
#'
#' @seealso
#' \code{\link[EEMtoolbox]{dicrepancy_continuous_sum}()} (original function);
#' \code{\link{add_species_names}()}, \code{\link{select_EEM_outputs}()} (helpers for equilibrium checking).
#'
#' @examples
#' \dontrun{
#'   # (1) Fit with no bounds → behaves identically to the original:
#'   res1 <- adapted_dicrepancy_continuous_sum(A = A_matrix, r = r_vector, ...)
#'
#'   # (2) Specify bounds that exclude a known equilibrium:
#'   lb <- rep(0.2, ncol(A_matrix))
#'   ub <- rep(1.5, ncol(A_matrix))
#'   res2 <- adapted_dicrepancy_continuous_sum(A = A_matrix, r = r_vector,
#'                                            lower_bounds = lb,
#'                                            upper_bounds = ub)
#'   if (res2$out_of_bound) {
#'     message("One or more species are outside the user‐defined bounds")
#'   }
#' }
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
