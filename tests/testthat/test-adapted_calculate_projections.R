# Test that projections starting at equilibrium remain constant
test_that("adapted_calculate_projections maintains equilibrium when initialized at steady state", {
  M <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  # Extract equilibrium
  eq_vals <- extract_eq(M)
  # Run projections with no interventions
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = eq_vals,
    t_window = c(0, 5),
    time_step_len = 0.1,
    model = "GLV",
    mode = "recruitment",
    multiplier = 1
  )
  n_steps <- length(seq(0, 5, by = 0.1))
  # All populations at each time should equal equilibrium
  for (i in seq_len(length(M[[1]]$growthrates))) {
    expect_true(all(abs(res[res$species == i & res$sim == 1,]$pop -
                          rep(eq_vals[[1]][i], length = n_steps)) < 1e-3))
    expect_true(all(abs(res[res$species == i & res$sim == 2,]$pop -
                          rep(eq_vals[[2]][i], length = n_steps)) < 1e-3))
  }
})

# Test events in adapted_calculate_projections
test_that("function triggers events at correct times", {
  M <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)

  #test recruitment events

#test init events triggered at correct times
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = rep(0,8),
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = "recruitment",
    init_intervention_amount = 5,
    init_intervention_timepoints = 0.5)
  expect_true(all(res[res$time == 0.5 & res$species == 1,]$pop == 0))
  expect_true(all(res[res$time == 0.6 & res$species == 1,]$pop > 0))
  #test sustain events triggered at correct times
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = rep(0,8),
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = "recruitment",
    sustain_intervention_amount = 5,
    sustain_intervention_timepoints = 0.5,
    sustain_intervention_threshold = 1000)
  expect_true(all(res[res$time == 0.5 & res$species == 1,]$pop == 0))
  expect_true(all(res[res$time == 0.6 & res$species == 1,]$pop > 0))
  #test threshold
  eq_vals <- extract_eq(M)
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = eq_vals,
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = "recruitment",
    sustain_intervention_amount = 5,
    sustain_intervention_timepoints = seq(0.5,1,0.1),
    sustain_intervention_threshold = min(eq_vals[[1]][1], eq_vals[[2]][1]) + 5)
  expect_true(all(round(res[res$time == 0.5 & res$species == 1,]$pop, 2) ==
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]), 2)))
  for (i in 1:2) {
    if (eq_vals[[i]][1] < min(eq_vals[[1]][1], eq_vals[[2]][1]) + 5) {
      expect_true(all(res[res$species == 1 & res$sim == i,]$pop <
                    min(eq_vals[[1]][1], eq_vals[[2]][1]) + 7))
    } else {
      expect_true(all(round(res[res$species == 1 & res$sim == i,]$pop, 2) ==
                    round(eq_vals[[i]][1], 2)))
    }
  }
#same with 2 species
  #test init events triggered at correct times
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = rep(0,8),
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = c("recruitment", "recruitment"),
    init_intervention_amount = c(5, 5),
    init_intervention_timepoints = list(0.5, 0.5),
    intro_species_index = c(1, 2),
    sustain_intervention_timepoints = list(NA, NA))
  expect_true(all(res[res$time == 0.5 & res$species == 1,]$pop == 0))
  expect_true(all(res[res$time == 0.6 & res$species == 1,]$pop > 0))
  expect_true(all(res[res$time == 0.5 & res$species == 2,]$pop == 0))
  expect_true(all(res[res$time == 0.6 & res$species == 2,]$pop > 0))

  #test removal events

  #test init events triggered at correct times
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = eq_vals,
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = "removal",
    init_intervention_amount = -5,
    init_intervention_timepoints = 0.5)
  expect_true(all(round(res[res$time == 0.5 & res$species == 1,]$pop, 2) ==
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]),2)))
  expect_true(all(round(res[res$time == 0.6 & res$species == 1,]$pop, 2) <
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]), 2)))
  #test sustain events triggered at correct times
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = eq_vals,
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = "removal",
    sustain_intervention_amount = -min(eq_vals[[1]][1], eq_vals[[2]][1])/2,
    sustain_intervention_timepoints = seq(0.5, 1, 0.1),
    sustain_intervention_threshold = 0)
  expect_true(all(round(res[res$time == 0.5 & res$species == 1,]$pop, 2) ==
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]), 2)))
  expect_true(all(round(res[res$time == 0.6 & res$species == 1,]$pop, 2) <
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]), 2)))
  expect_true(min(res[res$species == 1, ]$pop) == 0)
  #test threshold
  res <- adapted_calculate_projections(
    parameters = M,
    initial_condition = eq_vals,
    t_window = c(0, 1),
    time_step_len = 0.1,
    model = "GLV",
    scaled = FALSE,
    mode = "removal",
    sustain_intervention_amount = -min(eq_vals[[1]][1], eq_vals[[2]][1])/2,
    sustain_intervention_timepoints = seq(0.5, 1, 0.1),
    sustain_intervention_threshold = (min(eq_vals[[1]][1],
                                          eq_vals[[2]][1])/2) + 1)
  expect_true(all(round(res[res$time == 0.5 & res$species == 1,]$pop, 2) ==
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]),2)))
  expect_true(all(round(res[res$time == 0.6 & res$species == 1,]$pop, 2) <
                    round(c(eq_vals[[1]][1], eq_vals[[2]][1]), 2)))
  expect_true(all(res[res$species == 1, ]$pop > 0))
  })
