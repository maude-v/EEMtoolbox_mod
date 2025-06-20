#' @title extract_eq
#' @author Maude Vernet
#' @description
#' Extract equilibrium values from parameter sets of an EEM object
#' @param EEM A list of EEM objects
#' @param model The model used for the EEM objects. Default is "GLV".
#' @return A list of equilibrium values for each EEM object
#' @examples
#' out <- EEMtoolbox::EEM(matrix(c(-1, -1, 1, -1), ncol = 2), n_ensemble = 2)
#' equilibrium_values <- extract_eq(out)
#' @export
extract_eq <- function(EEM,
                       model = "GLV") {
  equilibrium_values <- vector(mode = "list", length = length(EEM))
  for (i in seq_len(length(EEM))) {
    if (model %in% c("GLV", "Gompertz")) {
      steady_state <- solve(EEM[[i]]$interaction_matrix,-EEM[[i]]$growthrates)
      if (model == "Gompertz") {
        steady_state <- exp(steady_state)
      }
    } else if (model == "Bimler-Baker") {
      r <- EEM[[i]]$growthrates
      A <- EEM[[i]]$interaction_matrix_alphas
      B <- EEM[[i]]$interaction_matrix_betas

      R <- r
      P <- diag(A)
      M <- A - diag(P)

      fn <- function(N) {
        output <- R*(1 - exp(-M %*% N - P)) + B %*% N #change to positive
        return(output)
      }
      sol <- nleqslv::nleqslv(rep(100,length(r)), fn) #we give a large positive warmstart
      steady_state <- sol$x
    }
    equilibrium_values[[i]] <- steady_state
  }
  return(equilibrium_values)
}
