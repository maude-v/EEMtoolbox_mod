# Test that all ensembles are selected when bounds are infinite (native mode)
test_that("All outputs selected for infinite bounds in native mode", {
  M <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  # Use infinite bounds to include all equilibria
  lower <- rep(-Inf, length(M[[1]]$growthrates))
  upper <- rep(Inf, length(M[[1]]$growthrates))
  outputs <- select_EEM_outputs(
    ensemble = M,
    target_lower = lower,
    target_upper = upper,
    mode = "native"
  )
  expect_equal(outputs$outputs_selected, seq_along(M))
  expect_length(outputs$outputs_unselected, 0)
})

# Test that no ensembles are selected when bounds exclude all equilibria (native mode)
test_that("No outputs selected for too-narrow bounds in native mode", {
  M <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  # Use very high lower bounds to exclude all equilibria
  lower <- rep(1e6, length(M[[1]]$growthrates))
  upper <- rep(1e7, length(M[[1]]$growthrates))
  # Capture printed message
  expect_output(
    {
      outputs <- select_EEM_outputs(
        ensemble = M,
        target_lower = lower,
        target_upper = upper,
        mode = "native"
      )
    },
    "No parameter sets found within the bounds."
  )
  expect_length(outputs$outputs_selected, 0)
  expect_equal(outputs$outputs_unselected, seq_along(M))
})

# Test disturbed mode: only the first species' equilibrium is checked
test_that("Disturbed mode selects ensembles based on first species bound", {
  M <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  # Determine reasonable infinite bounds for all but the first species
  lower <- c(-Inf, rep(-Inf, length(M[[1]]$growthrates) - 1))
  upper <- c(Inf, rep(Inf, length(M[[1]]$growthrates) - 1))
  outputs <- select_EEM_outputs(
    ensemble = M,
    target_lower = lower,
    target_upper = upper,
    mode = "disturbed",
    n_intro = 1
  )
  # Since bounds are infinite for species 1, all should be selected
  expect_equal(outputs$outputs_selected, seq_along(M))
  expect_length(outputs$outputs_unselected, 0)
})
