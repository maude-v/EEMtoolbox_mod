test_that("matrix gets rownames/colnames attached", {
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  new_names <- c("sp1", "sp2")
  outM <- add_species_names(output, new_names)[[1]]$interaction_matrix
  expect_equal(rownames(outM), new_names)
  expect_equal(colnames(outM), new_names)
})
