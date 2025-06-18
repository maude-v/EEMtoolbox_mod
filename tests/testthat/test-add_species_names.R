test_that("matrix gets rownames/colnames attached", {
  M <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  new_names <- c("sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7", "sp8")
  outM <- add_species_names(M, new_names)[[1]]$interaction_matrix
  expect_equal(rownames(outM), new_names)
  expect_equal(colnames(outM), new_names)
})
