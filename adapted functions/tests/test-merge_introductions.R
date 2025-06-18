# Ensure merge_introductions is available
test_that("Error when EEM_intros is a vector, not a pairlist", {
  M1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  M2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  flat_list <- list(M1, M2)
  expect_error(
    merge_introductions(flat_list),
    "parameter EEM_intros should be a pairlist"
  )
})

# Test mismatched projection lengths produce an error
test_that("Error when projections length differ", {
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 3)
  intros <- pairlist(E1, E2)
  expect_error(
    merge_introductions(intros),
    "same number of projections"
  )
})

# Test basic merge structure and dimensions
test_that("Merged output has correct structure and dimensions", {
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  intros <- pairlist(E1, E2)
  result <- merge_introductions(intros)

  # Should return a list of length equal to number of projections
  expect_length(result, length(E1))

  # Determine expected dimensions
  n_native <- ncol(E1[[1]]$interaction_matrix) - 1
  n_intro <- 2
  tot <- n_native + n_intro

  for (proj in result) {
    expect_true(is.list(proj))
    expect_true(all(c("growthrates", "interaction_matrix") %in% names(proj)))
    # growthrates length
    expect_length(proj$growthrates, tot)
    # interaction matrix dimensions
    mat <- proj$interaction_matrix
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), c(tot, tot))
  }
})

# Test default sign_interaction_intros = 0 replaces NAs with zeros
test_that("Default sign_interaction_intros replaces NA interactions with zero", {
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  intros <- pairlist(E1, E2)
  result <- merge_introductions(intros, sign_interaction_intros = 0)

  for (proj in result) {
    mat <- proj$interaction_matrix
    # No NA values remaining
    expect_false(any(is.na(mat)))
  }
})

# Test recycled vs updated modes produce consistent vs varying interaction terms
test_that("Recycled mode uses same interaction terms across projections", {
  set.seed(42)
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 3)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 3)
  intros <- pairlist(E1, E2)
  res_recycled <- merge_introductions(intros, sign_interaction_intros = 1, mode = "recycled")

  # For recycled, the interaction entry [1, 2] should be identical across projections
  term1 <- res_recycled[[1]]$interaction_matrix[1, 2]
  for (i in seq_along(res_recycled)) {
    expect_equal(res_recycled[[i]]$interaction_matrix[1, 2], term1)
  }
})

test_that("Updated mode uses new interaction terms each projection", {
  set.seed(42)
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 3)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 3)
  intros <- pairlist(E1, E2)
  res_updated <- merge_introductions(intros, sign_interaction_intros = 1, mode = "updated")

  # For updated, the same matrix entry should differ across projections
  vals <- sapply(res_updated, function(p) p$interaction_matrix[1, 2])
  expect_false(all(vals == vals[1]))
})

# Test that species names are preserved in dimnames
test_that("Species names are correctly set on merged interaction matrix", {
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  intros <- pairlist(E1, E2)
  result <- merge_introductions(intros)

  # Extract expected names
  n_native <- ncol(E1[[1]]$interaction_matrix) - 1
  intro_name1 <- colnames(E1[[1]]$interaction_matrix)[1]
  intro_name2 <- colnames(E2[[1]]$interaction_matrix)[1]
  native_names <- colnames(E1[[1]]$interaction_matrix)[2:(n_native + 1)]
  expected <- c(intro_name1, intro_name2, native_names)

  for (proj in result) {
    mat <- proj$interaction_matrix
    expect_equal(rownames(mat), expected)
    expect_equal(colnames(mat), expected)
  }
})

# Test that EEM1 and EEM2 have identical native columns
test_that("EEM1 and EEM2 share same native columns", {
  set.seed(123)
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- E1
  intros <- pairlist(E1, E2)
  # Compare original columns excluding first (intro)
  result <- merge_introductions(intros, sign_interaction_intros = 0)
  orig1 <- E1[[1]]$interaction_matrix[-1, -1]
  orig2 <- E2[[1]]$interaction_matrix[-1, -1]
  merged_mat <- result[[1]]$interaction_matrix[c(-1,-2),c(-1,-2)]
  expect_equal(as.vector(merged_mat), as.vector(orig1))
  expect_equal(as.vector(merged_mat), as.vector(orig2))
})

# Test mapping of first columns into merged output
test_that("Merged first column equals EEM1 first column", {
  set.seed(123)
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  intros <- pairlist(E1, E2)
  result <- merge_introductions(intros, sign_interaction_intros = 0)
  merged_mat <- result[[1]]$interaction_matrix
  orig1_mat <- E1[[1]]$interaction_matrix
  # merged_mat excluding row corresponding to EEM2
  comp <- merged_mat[-2, 1]
  expect_equal(comp, orig1_mat[,1])
})

# Test mapping of second column equals EEM2 first column
test_that("Merged second column equals EEM2 first column", {
  set.seed(123)
  E1 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  E2 <- EEMtoolbox::EEM(dingo_matrix, n_ensemble = 2)
  intros <- pairlist(E1, E2)
  result <- merge_introductions(intros, sign_interaction_intros = 0)
  merged_mat <- result[[1]]$interaction_matrix
  orig2_mat <- E2[[1]]$interaction_matrix
  comp <- merged_mat[-1, 2]
  expect_equal(comp, orig2_mat[,1])
})
