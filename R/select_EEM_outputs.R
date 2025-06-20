#' @title select_EEM_outputs
#' @author Maude Vernet
#' @description
#' This function filters a list of ecological ensemble model (EEM) outputs
#' based on whether the system has reached equilibrium and if the
#' equilibrium abundances of species fall within specified bounds.
#' @param ensemble an object of type `EEM`.
#' @param target_lower A numeric vector (or scalar) specifying the minimum
#'    acceptable equilibrium abundance for each species.
#' @param target_upper A numeric vector (or scalar) specifying the maximum
#'    acceptable equilibrium abundance for each species.
#' @param mode A character string indicating the mode of the simulation, either
#'   `"native"` or `"disturbed"`. This affects how the equilibrium is evaluated.
#' @param n_intro An integer specifying the number of introduced species in the
#'  simulation. This is used to label species in the output plot.
#'
#' @details
#' This function examines each element of `outputs` to (1) verify that the
#' system has reached equilibrium, and (2) check that every species’ equilibrium
#' abundance lies between `lower_bound` and `upper_bound`.
#'
#' @return A list containing two elements:
#'   - `outputs_selected`: Indices of the parameter sets that meet the criteria.
#'   - `outputs_unselected`: Indices of the parameter sets that do not meet the
#'   criteria.
#'
#' @examples
#' output <- EEMtoolbox::EEM(matrix(c(-1, -1, 1, -1), ncol = 2), n_ensemble = 2)
#' output <- select_EEM_outputs(output,
#'                              target_lower  = rep(0.1, 2),
#'                              target_upper  = rep(15, 2))

#' @export
select_EEM_outputs <- function(ensemble,
                               target_lower,
                               target_upper,
                               mode = "native",
                               n_intro = 1) {

  # Prepare an empty list to collect selected parameter sets
  selected <- list()
  unselected <- list()
  outputs <- list(outputs_selected = c(),
                  outputs_unselected = c())

  for (i in seq_along(ensemble)) {
    # Reconstruct the parameters: growth rates and interaction matrix.
    r <- ensemble[[i]]$growthrates
    A <- ensemble[[i]]$interaction_matrix

    # Compute the equilibrium: A %*% N + r = 0  =>  N = solve(A, -r)
    current_eq <- tryCatch(solve(A, -r),
                           error = function(e) rep(NA, length(r)))

    # Check for NA values (if the system is unsolvable, skip this set)
    if (any(is.na(current_eq))) {
      next
    }
    if (mode == "native") {
      # If the equilibrium is within the bounds, select this parameter set.
      if (length(which(current_eq > target_lower)) == length(r) &&
          length(which(current_eq < target_upper)) == length(r)) {
        selected[[length(selected) + 1]] <-
          list(ensemble[[i]],
               data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                          "species" = 1:length(r)))
        outputs$outputs_selected <-
          c(outputs$outputs_selected, i)
      } else {
        unselected[[length(unselected) + 1]] <-
          list(ensemble[[i]],
               data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                          "species" = 1:length(r)))
        outputs$outputs_unselected <-
          c(outputs$outputs_unselected, i)
      }
    }
    if (mode == "disturbed") {
      # If the equilibrium is within the bounds, select this parameter set.
      if (current_eq[1] > target_lower[1] &&
          current_eq[1] < target_upper[1]) {
        selected[[length(selected) + 1]] <-
          list(ensemble[[i]],
               data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                          "species" = 1:length(r)))
        outputs$outputs_selected <-
          c(outputs$outputs_selected, i)
      } else {
        unselected[[length(unselected) + 1]] <-
          list(ensemble[[i]],
               data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                          "species" = 1:length(r)))
        outputs$outputs_unselected <-
          c(outputs$outputs_unselected, i)
      }
    }
  }
  if (length(selected) == 0) {
    cat("No parameter sets found within the bounds.\n")
    df <- data.frame(equilibrium = numeric(), species = numeric())
    for (i in c(1:length(unselected))) {
      df <- data.frame(equilibrium = c(df$equilibrium,
                                       unselected[[i]][[2]]$equilibrium),
                       species = c(df$species,
                                   unselected[[i]][[2]]$species))
    }
  } else {
    cat("Selected", length(selected), "parameter sets.\n")
    df <- data.frame(equilibrium = numeric(), species = numeric())
    for (i in c(1:length(selected))) {
      df <- data.frame(equilibrium = c(df$equilibrium,
                                       selected[[i]][[2]]$equilibrium),
                       species = c(df$species,
                                   selected[[i]][[2]]$species))
    }
  }

  if (mode == "disturbed") {
    for (i in 1:length(df$species)) {
      if (df$species[i] %in% seq_len(n_intro)) {
        df$species[i] <- paste("introduced sp",
                               df$species[i], sep = " ")
      } else {
        df$species[i] <- paste("native sp",
                               as.numeric(df$species[i]) - n_intro, sep = " ")

      }
    }
  }

  uf <- data.frame(upper = target_upper, species = unique(df$species))
  lf <- data.frame(lower = target_lower, species = unique(df$species))

  a <- ggplot2::ggplot() +
    ggplot2::labs(y = "Computed equilibrium") +
    ggplot2::facet_wrap( ~ species, scales = "free") +
    ggplot2::geom_point(data = df,
                        ggplot2::aes(x = species, y = equilibrium,
                                     color = "equilibrium abundance")) +
    ggplot2::geom_point(data = uf,
                        ggplot2::aes(x = species, y = upper,
                                     color = "target upper")) +
    ggplot2::geom_point(data = lf,
                        ggplot2::aes(x = species, y = lower,
                                     color = "target lower")) +
    ggplot2::scale_color_manual(
      values = c("equilibrium abundance" = "forestgreen",
                 "target upper" = "red3",
                 "target lower" = "skyblue2")) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())


  print(a)
  return(outputs)
}
