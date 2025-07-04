#' @title adapted_calculate_projections
#' @author Maude Vernet
#' @description
#' Solves ODEs and plot solutions
#' @param parameters list like object of ensemble of parameters (outputs of  \link[EEMtoolbox]{EEM})
#' @param initial_condition vector of initial species abundances. If parameter scaled is TRUE, the parameter initial_condition should be scaled to steady state
#' @param t_window time window to solve ODE
#' @param time_step_len length of each time step, default = 0.01
#' @param model model representing species interactions. Default "GLV" (Generalized Lotka Volterra). options include "Bimler-Baker" and "Gompertz".
#' @param derivative derivative function. Default \link[EEMtoolbox]{derivative_func}
#' @param scaled Boolean indicating if projections should be scaled to steady state. If true, the parameter initial_condition should be scaled too. Default FALSE
#' @param species_names vector of strings for names of species. If NA plots only display species index number, . Default NA.
#' @param mode "recruitment" or "removal". Default "recruitment"
#' @param init_intervention_timepoints vector or list (if multiple introduced species) with the initial intervention timepoints for each introduced species, default = NA
#' @param init_intervention_amount vector with the number of individuals in initial interventions for each introduced species, default = 0
#' @param sustain_intervention_amount vector with the number of individuals in subsequent interventions for each introduced species, default = 0
#' @param sustain_intervention_timepoints vector or list (if multiple introduced species) with the subsequent intervention timepoints for each introduced species, default = NA
#' @param sustain_intervention_threshold vector or list (if multiple introduced species) with the abundance threshold where subsequent interventions stop for each introduced species, default = NA
#' @param intro_species_index vector indicating which species in the index is/are the introduced one(s), default = 1
#' @param multiplier multiplier for the initial condition, default = 1
#' @param extinction_threshold threshold below which species are considered extinct, default = 0
#' @examples
#' library(EEMtoolbox)
#' output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol =2), n_ensemble = 2)
#' adapted_calculate_projections(output,
#'                               initial_condition = rep(1,2),
#'                               t_window = c(0,1),
#'                               mode = "recruitment",
#'                               init_intervention_amount = 2,
#'                               init_intervention_timepoints = c(1, 2))
#'
#' @return dataset of species abundances over time
#' @export
adapted_calculate_projections <-
  function(parameters,
           initial_condition,
           t_window,
           time_step_len=0.01,
           model = "GLV",
           derivative=EEMtoolbox::derivative_func,
           scaled=FALSE,
           species_names=NA,
           # Recruitment parameters:
           mode = "recruitment",
           init_intervention_amount = 0,
           init_intervention_timepoints = NA,
           sustain_intervention_amount = 0,
           sustain_intervention_timepoints = NA,
           sustain_intervention_threshold = NA,
           intro_species_index = 1,
           multiplier = 1,
           extinction_threshold = 0) {

    if (length(intro_species_index) == 1) {
      init_intervention_timepoints <- list(init_intervention_timepoints)
      sustain_intervention_timepoints <- list(sustain_intervention_timepoints)
    }

    # We'll pass the recruitment schedule via the parameters list.
    # Define the intervention times over the time window:
    recruitment_times <-
      round(seq(from = t_window[1],
                to = t_window[2],
                by = time_step_len),
            nchar(strsplit(
              as.character(time_step_len), "\\.")[[1]][2]))

    recruitment_event <- function(time, y, pars) {
      # Check if current time is an initial intervention time:
      # If so, add the initial intervention amount to the introduced species:
      for (i in seq_along(pars$intro_species_index)) {
        if (time %in% pars$init_intervention_timepoints[[i]]) {
          y[pars$intro_species_index[i]] <-
            #in this case, y is the current vector of species abundance at a given time
            y[pars$intro_species_index[i]] + pars$init_intervention_amount[i]
        }
        if (mode[i] == "removal") {
          if (time %in% pars$sustain_intervention_timepoints[[i]] &&
            y[pars$intro_species_index[i]] > pars$sustain_intervention_threshold[i]) {
            y[pars$intro_species_index[i]] <-
              y[pars$intro_species_index[i]] + pars$sustain_intervention_amount[i]
          }
        } else if (mode[i] == "recruitment") {
          # Check if current time is a sustaining intervention time AND abundance is below threshold:
          # If so, add the sustaining intervention amount to the introduced species:
          if (time %in% pars$sustain_intervention_timepoints[[i]] &&
            y[pars$intro_species_index[i]] < pars$sustain_intervention_threshold[i]) {
            y[pars$intro_species_index[i]] <-
              y[pars$intro_species_index[i]] + pars$sustain_intervention_amount[i]
          }
        }
      }
      # extinction
      y[y <= pars$extinction_threshold] <- 0
      return(y)
    }

    # Create a recruitment parameter list to pass to the ODE solver. This way, we
    # don't need to write them one after the other but rather keep them all
    # together.
    recruitment_pars <-
      list(init_intervention_timepoints = init_intervention_timepoints,
           sustain_intervention_timepoints = sustain_intervention_timepoints,
           init_intervention_amount = init_intervention_amount,
           sustain_intervention_amount = sustain_intervention_amount,
           sustain_intervention_threshold = sustain_intervention_threshold,
           intro_species_index = intro_species_index)

    ode_solve_it <- function(pars,
                             initial_condition,
                             t_window,
                             time_step_len,
                             model,
                             derivative,
                             scaled,
                             recruitment_event, #added -> function defined earlier
                             recruitment_times, #added -> we defined earlier. Do we need this?
                             recruitment_pars, #added -> all the parameters related to recruitement, packed together.
                             extinction_threshold) {
      #if scaled: find steady state, initial condition is scaled to steady state
      if (scaled) { #if scaled == TRUE
        if (model %in% c("GLV", "Gompertz")) {
          steady_state <- solve(pars$interaction_matrix,-pars$growthrates)
          if (model == "Gompertz") {
            steady_state <- exp(steady_state)
          }
        } else if (model == "Bimler-Baker") {
          r <- pars$growthrates
          A <- pars$interaction_matrix_alphas
          B <- pars$interaction_matrix_betas

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
        initial <- initial_condition * steady_state #initial condition is scaled to steady state -> interpret initial as a multiplier
      } else {
        initial <- initial_condition #non-scaled to steady state
      }

      ## solve ODE
      if (model %in% c("GLV", "Gompertz")) {
        projections <-
          as.data.frame(adapted_ode_solve( #use adapted ode_solve function
            interaction_matrix_value = pars$interaction_matrix,
            growth_rate = pars$growthrates,
            initial,
            t_window,
            time_step_len,
            model,
            derivative,
            recruitment_pars = recruitment_pars,
            recruitment_event = recruitment_event,
            recruitment_times = recruitment_times,
            extinction_threshold = extinction_threshold)) #added
      } else if (model == "Bimler-Baker") {
        projections <-
          as.data.frame(adapted_ode_solve( #use adapted ode_solve function
            interaction_matrix_value = list(pars$interaction_matrix_alphas,
                                            pars$interaction_matrix_betas),
            growth_rate = pars$growthrates,
            initial,
            t_window,
            time_step_len,
            model,
            derivative,
            recruitment_pars = recruitment_pars,
            recruitment_event = recruitment_event,
            recruitment_times = recruitment_times,
            extinction_threshold = extinction_threshold)) #added
      }

      if (scaled) { #if scaled to steady state, scale abundances
        projections[,-1] <- projections[,-1]/matrix(steady_state,
                                                    ncol = length(steady_state),
                                                    nrow = nrow(projections),
                                                    byrow = TRUE)
      }
      return(projections)
    }

    if (isFALSE(is.list(initial_condition))) {
      list_init <- vector(mode = "list", length = length(parameters))
      for (i in seq_len(length(parameters))) {
        list_init[[i]] <- initial_condition
      }
      initial_condition <- list_init
    }

    abundance <- lapply(seq_len(length(parameters)),
                        function(i) {
                          cat("\r","Solving ODE for parameter set", i, "/", length(parameters))
                          ode_solve_it(parameters[[i]],
                                       model = model,
                                       initial_condition = initial_condition[[i]],
                                       t_window = t_window,
                                       time_step_len = time_step_len,
                                       derivative,
                                       scaled,
                                       recruitment_event = recruitment_event, # added
                                       recruitment_times = recruitment_times, # added
                                       recruitment_pars = recruitment_pars,
                                       extinction_threshold = extinction_threshold) # added
                        })

    abundance <- dplyr::bind_rows(abundance, .id = "sim")

    if (!is.na(species_names[1])) {
      names(abundance) <- c("sim", "time", species_names)
    }

    ## remove projections that could not be solved ####
    remove_indexes <-  dplyr::group_by(abundance, sim)
    remove_indexes <-  dplyr::summarise(remove_indexes, max_time = max(time))
    remove_indexes <-  dplyr::filter(remove_indexes, max_time < t_window[2])
    remove_indexes <- remove_indexes$sim

    if (length(remove_indexes) > 0) {
      print("Following parameter sets could not be solved and were removed:")
      print(remove_indexes)

      #remove parameter sets#
      # abundance <- dplyr::filter(abundance, !(sim %in% remove_indexes)) #this line was already in green
      abundance <- abundance[which(!(abundance$sim %in% remove_indexes)),]
    }

    #pivot for plotting
    abundance <- tidyr::pivot_longer(abundance,
                                     !c(time,sim), #select all columns except time and sim
                                     names_to = c("species"),
                                     values_to = "pop")
    if (scaled == FALSE) {
      abundance$pop <- abundance$pop*multiplier #multiply by the multiplier, dropped the rest because we replaced it with the root
    }

    return(abundance)
  }
