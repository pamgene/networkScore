make_results_csv <- function(dir, comparisons, score_sig_network, score_sig_network_inv) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(
    data.frame(
      Comparison = comparisons, score_sig_network = score_sig_network,
      score_sig_network_inv = score_sig_network_inv
    ),
    file.path(dir, "results.csv")
  )
}

test_that("plot_score_comparison facets by alias, in the order given, with all comparisons present", {
  dir_a <- file.path(tempdir(), "compare_a")
  dir_b <- file.path(tempdir(), "compare_b")
  unlink(c(dir_a, dir_b), recursive = TRUE)
  make_results_csv(dir_a, c("condA", "condB"), c(0.1, 0.2), c(0.3, 0.4))
  make_results_csv(dir_b, c("condA", "condB"), c(0.5, 0.6), c(0.7, 0.8))

  p <- plot_score_comparison(c(baseline = dir_a, with_KL = dir_b))

  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$alias), c("baseline", "with_KL"))
  expect_setequal(unique(p$data$Comparison), c("condA", "condB"))
  expect_setequal(unique(p$data$metric), c("score_sig_network", "score_sig_network_inv"))
  # 2 folders x 2 comparisons x 2 metrics = 8 long-format rows
  expect_equal(nrow(p$data), 8)
})

test_that("plot_score_comparison finds results.csv recursively under a multi-combination folder", {
  base_dir <- file.path(tempdir(), "compare_nested")
  unlink(base_dir, recursive = TRUE)
  make_results_csv(file.path(base_dir, "b50"), "condA", 0.1, 0.2)
  make_results_csv(file.path(base_dir, "b75"), "condA", 0.3, 0.4)

  p <- plot_score_comparison(c(run1 = base_dir))

  # Both b50 and b75 rows found and combined -- 2 rows x 2 metrics = 4
  expect_equal(nrow(p$data), 4)
})

test_that("plot_score_comparison drops comparisons not common to every folder, with a message", {
  dir_a <- file.path(tempdir(), "compare_partial_a")
  dir_b <- file.path(tempdir(), "compare_partial_b")
  unlink(c(dir_a, dir_b), recursive = TRUE)
  make_results_csv(dir_a, c("condA", "condB"), c(0.1, 0.2), c(0.3, 0.4))
  make_results_csv(dir_b, c("condA", "condC"), c(0.5, 0.6), c(0.7, 0.8))

  expect_message(
    p <- plot_score_comparison(c(x = dir_a, y = dir_b)),
    "condB|condC"
  )
  expect_setequal(unique(p$data$Comparison), "condA")
})

test_that("plot_score_comparison errors clearly when no comparison is shared across folders", {
  dir_a <- file.path(tempdir(), "compare_disjoint_a")
  dir_b <- file.path(tempdir(), "compare_disjoint_b")
  unlink(c(dir_a, dir_b), recursive = TRUE)
  make_results_csv(dir_a, "condA", 0.1, 0.2)
  make_results_csv(dir_b, "condZ", 0.5, 0.6)

  expect_error(
    plot_score_comparison(c(x = dir_a, y = dir_b)),
    "No Comparison is common"
  )
})

test_that("plot_score_comparison errors clearly when a folder has no results.csv", {
  empty_dir <- file.path(tempdir(), "compare_empty")
  unlink(empty_dir, recursive = TRUE)
  dir.create(empty_dir)

  expect_error(
    plot_score_comparison(c(x = empty_dir)),
    "No results.csv found"
  )
})

test_that("plot_score_comparison requires named respaths", {
  expect_error(
    plot_score_comparison(c(tempdir())),
    "named"
  )
})

test_that("plot_score_comparison errors clearly when a requested metric column is missing", {
  dir_a <- file.path(tempdir(), "compare_missing_metric")
  unlink(dir_a, recursive = TRUE)
  dir.create(dir_a)
  readr::write_csv(data.frame(Comparison = "condA", score_sig_network = 0.1), file.path(dir_a, "results.csv"))

  expect_error(
    plot_score_comparison(c(x = dir_a)),
    "score_sig_network_inv"
  )
})
