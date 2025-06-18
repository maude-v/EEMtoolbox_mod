# Sample projections data for a single species at two time points
proj_df <- data.frame(
  sim = c(1, 1, 1, 1, 1, 1, 1, 1,
          2, 2, 2, 2, 2, 2, 2, 2),
  time = c(0, 0, 1, 1, 2, 2, 3, 3,
           0, 0, 1, 1, 2, 2, 3, 3),
  species = rep(c("S", "M"), 8),
  pop = c(1, 1, 3, 3, 2, 2, 4, 4,
          2, 2, 6, 6, 4, 4, 12, 12))

# Test aggregation and multiplier
test_that("adapted_plot_projections correctly aggregates pop and applies multiplier", {
  p <- adapted_plot_projections(projections = proj_df,
                                title = NA,
                                multiplier = 3,
                                scaled = TRUE)
  # Check ggplot object
  expect_s3_class(p, "ggplot")

  # Build plot data
  bd <- ggplot2::ggplot_build(p)
  # Ribbon is first layer
  ribbon <- bd$data[[1]]
  # Line is second layer
  line <- bd$data[[2]]

  # Expect two x values: 0, 1, 2, 3
  expect_equal(unique(ribbon$x), c(0, 1, 2, 3))
  # Aggregation: time=0 pop values 1,3 -> median=2, min=1, max=3
  # After multiplier=3: median=6, ymin=3, ymax=9
  # time=1 pop values 2,4 -> median=3, min=2, max=4
  # After multiplier: median=9, ymin=6, ymax=12
  # ribbon$ymin corresponds to lower bound (min * multiplier)
  expect_equal(ribbon$ymin, rep(c(3, 9, 6, 12), 2))
  # ribbon$ymax corresponds to upper bound (max * multiplier)
  expect_equal(ribbon$ymax, rep(c(6, 18, 12, 36), 2))

  # line$y corresponds to median * multiplier
  expect_equal(line$y, rep(c(4.5, 13.5, 9, 24), 2))
})

# Test title parameter is applied
test_that("adapted_plot_projections sets plot title when provided", {
  p <- adapted_plot_projections(
    projections = proj_df,
    title = "My Title",
    multiplier = 1,
    scaled = TRUE
  )
  expect_equal(p$labels$title, "My Title")
})

# Test scaled parameter controls facet scales
test_that("adapted_plot_projections uses free scales when scaled=TRUE and fixed when scaled=FALSE", {
  p_free <- adapted_plot_projections(projections = proj_df,
    title = NA,
    multiplier = 1,
    scaled = TRUE)
  expect_equal(p_free$facet$params$free$x, TRUE)
  expect_equal(p_free$facet$params$free$y, TRUE)

  p_fixed <- adapted_plot_projections(projections = proj_df,
    title = NA,
    multiplier = 1,
    scaled = FALSE)
  expect_equal(p_fixed$facet$params$free$x, FALSE)
  expect_equal(p_fixed$facet$params$free$y, FALSE)
})
