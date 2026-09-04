test_that("plot_network_metric_hist_facet writes one PNG per metric x cutoff, faceted by condition", {
  respath <- file.path(tempdir(), "plot_facet_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  metrics_df <- data.frame(
    condition = rep(c("A", "B"), each = 20),
    perc_cutoff = 0.7,
    metric = "rel_med_path",
    value = c(rnorm(20, 0.5), rnorm(20, 0.6))
  )
  obs_df <- data.frame(condition = c("A", "B"), perc_cutoff = 0.7, metric = "rel_med_path", value = c(0.1, 0.9))

  plot_network_metric_hist_facet(metrics_df, obs_df, respath)

  expect_true(file.exists(file.path(respath, "Network_metrics_plots", "rel_med_path_cutoff0.7.png")))
})

test_that("plot_network_metric_hist_facet handles an empty metrics data frame without erroring", {
  respath <- file.path(tempdir(), "plot_facet_empty_test")
  dir.create(respath, showWarnings = FALSE)
  expect_no_error(plot_network_metric_hist_facet(data.frame(), data.frame(), respath))
})
