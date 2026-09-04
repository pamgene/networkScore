test_that("initialize_or_read_df starts empty, then picks up a written file on resume", {
  respath <- file.path(tempdir(), "init_read_test")
  dir.create(respath, showWarnings = FALSE)
  unlink(file.path(respath, "results.csv"))

  fresh <- initialize_or_read_df(respath, "results")
  expect_equal(nrow(fresh), 0)

  readr::write_csv(data.frame(a = 1:2), file.path(respath, "results.csv"))
  resumed <- initialize_or_read_df(respath, "results")
  expect_equal(nrow(resumed), 2)
})

test_that("finalize_golden_score_results writes csvs, deletes temp files, and returns results+logs", {
  respath <- file.path(tempdir(), "finalize_test")
  dir.create(respath, showWarnings = FALSE)

  temp_file <- file.path(respath, "temp_network_condA_0.7.txt")
  writeLines("0.05", temp_file)

  out <- finalize_golden_score_results(
    results = data.frame(condition = "condA", score_sig_network = 0.05),
    all_metrics_df = data.frame(),
    all_obs_metrics_df = data.frame(),
    temp_files = temp_file,
    respath = respath,
    loop_start = Sys.time(),
    logs = character()
  )

  expect_true(file.exists(file.path(respath, "results.csv")))
  expect_false(file.exists(temp_file))
  expect_equal(out$results$condition, "condA")
  expect_true(length(out$logs) > 0)
})

test_that("reconstruct_golden_score_from_temp and find_missing_golden_score_combinations agree with what's on disk", {
  respath <- file.path(tempdir(), "reconstruct_test")
  dir.create(respath, showWarnings = FALSE)
  writeLines("0.02", file.path(respath, "temp_network_condA_0.7.txt"))

  reconstructed <- reconstruct_golden_score_from_temp(respath, conditions = c("condA", "condB"), perc_cutoffs = 0.7)
  expect_equal(reconstructed$score_sig_network[reconstructed$condition == "condA"], 0.02)
  expect_true(is.na(reconstructed$score_sig_network[reconstructed$condition == "condB"]))

  missing <- find_missing_golden_score_combinations(respath, conditions = c("condA", "condB"), perc_cutoffs = 0.7)
  expect_equal(missing$condition, "condB")
})
