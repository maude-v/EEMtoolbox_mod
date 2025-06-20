# Test that native parameters are preserved and introduced parameters are within bounds
test_that("add_introduced_species preserves native parameters and applies bounds and signs correctly", {
  set.seed(123)
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  n_native <- length(output[[1]]$growthrates)
  lower <- 1
  upper <- 5
  self_sign <- -1
  row_signs <- rep(1, n_native)
  col_signs <- rep(-1, n_native)
  k <- 2

  ext <- add_introduced_species(
    native_parameters = output,
    introduced_lower_bound_growth_rate = lower,
    introduced_upper_bound_growth_rate = upper,
    introduced_self_sign = self_sign,
    introduced_row_signs = row_signs,
    introduced_col_signs = col_signs,
    introduced_k = k
  )

  # Should return same number of parameter sets
  expect_length(ext, length(output))

  for (i in seq_along(ext)) {
    e <- ext[[i]]
    # Check growthrates
    expect_length(e$growthrates, n_native + 1)
    # Native growthrates unchanged
    expect_equal(e$growthrates[-1], output[[i]]$growthrates)
    # Introduced growthrate within bounds
    expect_true(e$growthrates[1] >= lower && e$growthrates[1] <= upper)

    # Check interaction matrix dimensions
    mat <- e$interaction_matrix
    expect_equal(dim(mat), c(n_native + 1, n_native + 1))

    # Native-native block unchanged
    expect_equal(
      as.vector(mat[2:(n_native + 1), 2:(n_native + 1)]),
      as.vector(output[[i]]$interaction_matrix)
    )

    # Check row and column signs for interactions
    expect_equal(as.vector(sign(mat[1, 2:(n_native + 1)])), row_signs)
    expect_equal(as.vector(sign(mat[2:(n_native + 1), 1])), col_signs)

    # No NA values remain
    expect_false(any(is.na(mat)))
  }
})

# Test that the equilibrium of the introduced species equals k
test_that("introduced species equilibrium equals specified k", {
  set.seed(456)
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  n_native <- length(output[[1]]$growthrates)
  lower <- 0.5
  upper <- 1.5
  self_sign <- 1
  row_signs <- rep(0, n_native)
  col_signs <- rep(0, n_native)
  k <- 3

  ext <- add_introduced_species(native_parameters = output,
                                introduced_lower_bound_growth_rate = lower,
                                introduced_upper_bound_growth_rate = upper,
                                introduced_self_sign = self_sign,
                                introduced_row_signs = row_signs,
                                introduced_col_signs = col_signs,
                                introduced_k = k)

  eqs <- extract_eq(ext, model = "GLV")
  for (i in seq_along(eqs)) {
    expect_equal(as.vector(eqs[[i]][1]), k, tolerance = 1e-6)
  }
})
