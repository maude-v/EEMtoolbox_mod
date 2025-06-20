# Ensure merge_introductions is available
test_that("Error when EEM_intros is a vector, not a pairlist", {
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  flat_list <- list(output1, output2)
  expect_error(
    merge_introductions(flat_list),
    "parameter EEM_intros should be a pairlist"
  )
})

# Test mismatched projection lengths produce an error
test_that("Error when projections length differ", {
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 3)
  intros <- pairlist(output1,output2)
  expect_error(
    merge_introductions(intros),
    "same number of projections"
  )
})

# Test basic merge structure and dimensions
test_that("Merged output has correct structure and dimensions", {
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  result <- merge_introductions(intros)

  # Should return a list of length equal to number of projections
  expect_length(result, length(output1))

  # Determine expected dimensions
  n_native <- ncol(output1[[1]]$interaction_matrix) - 1
  n_intro <- 2
  tot <- n_native + n_intro

  for (i in seq_along(result)) {
    expect_true(is.list(result[[i]]))
    expect_true(all(c("growthrates", "interaction_matrix") %in%
                      names(result[[i]])))
    # growthrates length
    expect_length(result[[i]]$growthrates, tot)
    # interaction matrix dimensions
    mat <- result[[i]]$interaction_matrix
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), c(tot, tot))
  }
})

# Test default sign_interaction_intros = 0 replaces NAs with zeros
test_that("Default sign_interaction_intros replaces NA interactions with zero", {
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  result <- merge_introductions(intros, sign_interaction_intros = 0)

  for (i in seq_along(result)) {
    mat <- result[[i]]$interaction_matrix
    # No NA values remaining
    expect_false(any(is.na(mat)))
  }
})

# Test recycled vs updated modes produce consistent vs varying interaction terms
test_that("Recycled mode uses same interaction terms across projections", {
  set.seed(42)
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  res_recycled <- merge_introductions(intros,
                                      sign_interaction_intros = 1,
                                      mode = "recycled")

  # For recycled, the interaction entry [1, 2] should be identical across projections
  term1 <- res_recycled[[1]]$interaction_matrix[1, 2]
  for (i in seq_along(res_recycled)) {
    expect_equal(res_recycled[[i]]$interaction_matrix[1, 2], term1)
  }
})

test_that("Updated mode uses new interaction terms each projection", {
  set.seed(42)
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  res_updated <- merge_introductions(intros,
                                     sign_interaction_intros = 1,
                                     mode = "updated")

  # For updated, the same matrix entry should differ across projections
  vals <- sapply(res_updated, function(p) p$interaction_matrix[1, 2])
  expect_false(all(vals == vals[1]))
})

# Test that species names are preserved in dimnames
test_that("Species names are correctly set on merged interaction matrix", {
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output1 <- add_species_names(output1, c("Native1", "Native2"))
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- add_species_names(output1, c("Native3", "Native4"))
  intros <- pairlist(output1, output2)
  result <- merge_introductions(intros)

  # Extract expected names
  n_native <- ncol(output1[[1]]$interaction_matrix) - 1
  intro_name1 <- colnames(output1[[1]]$interaction_matrix)[1]
  intro_name2 <- colnames(output2[[1]]$interaction_matrix)[1]
  native_names <- colnames(output1[[1]]$interaction_matrix)[2:(n_native + 1)]
  expected <- c(intro_name1, intro_name2, native_names)

  for (i in seq_along(result)) {
    mat <- result[[i]]$interaction_matrix
    expect_equal(rownames(mat), expected)
    expect_equal(colnames(mat), expected)
  }
})

# Test that EEM1 and EEM2 have identical native columns
test_that("EEM1 and EEM2 share same native columns", {
  set.seed(123)
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  # Compare original columns excluding first (intro)
  result <- merge_introductions(intros, sign_interaction_intros = 0)
  orig1 <- output1[[1]]$interaction_matrix[-1, -1]
  orig2 <- output2[[1]]$interaction_matrix[-1, -1]
  merged_mat <- result[[1]]$interaction_matrix[c(-1,-2),c(-1,-2)]
  expect_equal(as.vector(merged_mat), as.vector(orig1))
})

# Test mapping of first columns into merged output
test_that("Merged first column equals EEM1 first column", {
  set.seed(123)
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  result <- merge_introductions(intros, sign_interaction_intros = 0)
  merged_mat <- result[[1]]$interaction_matrix
  orig1_mat <- output1[[1]]$interaction_matrix
  # merged_mat excluding row corresponding to EEM2
  comp <- merged_mat[-2, 1]
  expect_equal(comp, orig1_mat[,1])
})

# Test mapping of second column equals EEM2 first column
test_that("Merged second column equals EEM2 first column", {
  set.seed(123)
  output1 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  output2 <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol = 2), n_ensemble = 2)
  intros <- pairlist(output1, output2)
  result <- merge_introductions(intros, sign_interaction_intros = 0)
  merged_mat <- result[[1]]$interaction_matrix
  orig2_mat <- output2[[1]]$interaction_matrix
  comp <- merged_mat[-1, 2]
  expect_equal(comp, orig2_mat[,1])
})
