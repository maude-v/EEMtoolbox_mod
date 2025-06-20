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
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2),
                            n_ensemble = 2,
                            disc_func = function(data) {
                              adapted_discrepancy_continuous_sum(
                                data,
                                target_lower = rep(1, 2),
                                target_upper = rep(15, 2))})
  eq <- extract_eq(output)
  target_lower <- rep(1, 2)
  target_upper <- rep(15, 2)
  for (i in seq_along(output)) {
  jacobian <- output[[i]]$interaction_matrix*matrix(eq[[i]],
                                                    ncol = 2,
                                                    nrow = 2)
  diag(jacobian) <- diag(jacobian) +
    output[[i]]$interaction_matrix %*% eq[[i]] +
    output[[i]]$growthrates
  stability_eigenvalues <- Re(eigen(jacobian)$values)
  res <- EEMtoolbox::adapted_discrepancy_continuous_sum(
    data = c(eq[[i]],
             stability_eigenvalues),
    target_lower = target_lower,
    target_upper = target_upper)
  expect_true(res == 0)
  }
})
