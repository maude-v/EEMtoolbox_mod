test_that("test discrepancy function", {

  test_vect <- rep(0,8)
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect), 0)

  test_vect2 <- c(rep(1,4), rep(-1,4)) #test stability and feasibility
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect2), 0)

  test_vect3 <- c(rep(-1,4), rep(-1,4)) #infeasible
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect3), 4)

  test_vect4 <- c(rep(1,4), rep(1,4)) #unstable
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect4), 4)

  test_vect5 <- c(rep(-1,4), rep(1,4)) #infeasible and unstable
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect5), 8)
})

test_that("bounds enforcement", {
  A_toy <- matrix(c(2, 0,
                    0, 2), nrow = 2)
  r_toy <- c(1, 1)
  # Solve for equilibrium: x* = A^{-1} r = 0.5, 0.5
  # Now impose a lower bound > 0.5, e.g. lower_bounds = c(0.7, 0.7)
  lower <- c(0.7, 0.7)
  upper <- NULL
  res <-
    EEMtoolbox::adapted_dicrepancy_continuous_sum(A = A_toy,
                                                  r = r_toy,
                                                  lower_bounds = lower,
                                                  upper_bounds = upper,
                                                  stop_on_out_of_bound = FALSE)
  expect_true(res$out_of_bound)
  # Both species violate the lower bound:
  expect_setequal(res$bounds_violation_indices, c(1L, 2L))
})
