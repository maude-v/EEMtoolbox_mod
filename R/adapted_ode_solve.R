#' @title adapted_ode_solve
#' @author Maude Vernet
#' @description
#' This function is a wrapper for the `deSolve::ode` function, which solves
#' ordinary differential equations (ODEs) for ecological models.
#' It allows for the inclusion of recruitment events and extinction thresholds,
#' making it suitable for ecological modeling scenarios.
#' It is adapted from the `EEMtoolbox::ode_solve` function.
#' @param initial_condition A numeric vector representing the initial conditions
#' of the system.
#' @param interaction_matrix_value A matrix representing the interaction
#' coefficients between species.
#' @param growth_rate A numeric vector representing the growth rates of the
#' species.
#' @param t_window A numeric vector of length 2 representing the time window for
#' the simulation (start and end times).
#' @param time_step_len A numeric value representing the length of each time
#' step in the simulation.
#' @param model A character string specifying the model type (default is "GLV"
#' for Generalized Lotka-Volterra).
#' @param derivative A function that computes the derivatives of the system
#' (default is `EEMtoolbox::derivative_func`).
#' @param recruitment_pars A list of parameters related to recruitment events.
#' @param recruitment_event A function that defines the recruitment event.
#' @param recruitment_times A numeric vector specifying the times at which
#' recruitment events occur.
#' @param extinction_threshold A numeric value representing the threshold below
#' which species are considered extinct.
#' @return A data frame containing the results of the ODE solution, including
#' time steps and species populations.
#' @author Maude Vernet
#' @export
adapted_ode_solve <- function(initial_condition,
                              interaction_matrix_value,
                              growth_rate,
                              t_window,
                              time_step_len,
                              model = "GLV",
                              derivative = EEMtoolbox::derivative_func,
                              recruitment_pars,
                              recruitment_event,
                              recruitment_times,
                              extinction_threshold){

  time_steps <- round(seq(from = t_window[1],
                          to = t_window[2],
                          by = time_step_len),
                      nchar(strsplit(
                        as.character(time_step_len), "\\.")[[1]][2])) #vector of time steps -> to make sure it is actually the correct sequence, there is a bug with the seq() of r

  pars  <- c(list(interaction_matrix_value = interaction_matrix_value,
                  growth_rate  = growth_rate,
                  model = model),
             recruitment_pars, #added, parameters related to recruitment
             list(extinction_threshold = extinction_threshold))

  out <- deSolve::ode(y = initial_condition,
                      times = time_steps,
                      func = derivative,
                      parms = pars,
                      events = list(func = recruitment_event,
                                    time = recruitment_times)) # added

  return(out)
}
