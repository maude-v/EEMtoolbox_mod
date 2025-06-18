# Test GLV equilibrium calculation
test_that("extract_eq computes GLV equilibrium correctly", {
  r <- c(-1, -4)
  A <- diag(c(1, 2))
  EEM <- list(list(growthrates = r, interaction_matrix = A))
  eq_vals <- extract_eq(EEM)
  expect_equal(eq_vals[[1]], c(1, 2))
})

# Test Gompertz equilibrium calculation
test_that("extract_eq computes Gompertz equilibrium correctly", {
  r <- c(-1, -4)
  A <- diag(c(1, 2))
  EEM <- list(list(growthrates = r, interaction_matrix = A))
  eq_vals <- extract_eq(EEM, model = "Gompertz")
  expect_equal(eq_vals[[1]], exp(c(1, 2)))
})

# Test Bimler-Baker returns correct length and type
test_that("extract_eq returns numeric vector of correct length for Bimler-Baker model", {
  r <- c(0.5, 1.0)
  A <- diag(c(1, 1))
  B <- matrix(0, nrow = 2, ncol = 2)
  # Create EEM object with alphas and betas
  EEM_BB <- list(list(
    growthrates = r,
    interaction_matrix_alphas = A,
    interaction_matrix_betas = B
  ))
  eq_vals <- extract_eq(EEM_BB, model = "Bimler-Baker")
  expect_true(is.numeric(eq_vals[[1]]))
  expect_length(eq_vals[[1]], length(r))
})

# Test default model is GLV
test_that("extract_eq uses GLV by default", {
  r <- c(-2, -3)
  A <- diag(c(2, 3))
  EEM <- list(list(growthrates = r, interaction_matrix = A))
  eq_default <- extract_eq(EEM)
  eq_glv <- extract_eq(EEM, model = "GLV")
  expect_equal(eq_default, eq_glv)
})
