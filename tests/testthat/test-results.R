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

test_that("drop_constant_grid_columns drops only the grid columns that don't vary", {
  results <- data.frame(
    Comparison = c("A", "B", "A", "B"),
    spec_cutoff = c(1, 1, 1, 1), # constant -- dropped
    perc_cutoff = c(0.5, 0.5, 0.7, 0.7), # varies -- kept
    b = c(50, 50, 75, 75), # varies -- kept
    rank_uka_abs = c(TRUE, TRUE, TRUE, TRUE), # constant -- dropped
    ppi_network = c("v12", "v12", "v12", "v12"), # constant -- dropped
    score_sig_network = c(0.1, 0.2, 0.3, 0.4) # not a grid column -- always kept
  )

  out <- drop_constant_grid_columns(results)

  expect_equal(colnames(out), c("Comparison", "perc_cutoff", "b", "score_sig_network"))
  expect_equal(nrow(out), 4)
})

test_that("drop_constant_grid_columns keeps every grid column when everything varies", {
  results <- data.frame(spec_cutoff = c(0.5, 0.7), perc_cutoff = c(0.5, 0.7), b = c(1, 2))
  out <- drop_constant_grid_columns(results)
  expect_equal(colnames(out), c("spec_cutoff", "perc_cutoff", "b"))
})

test_that("drop_constant_grid_columns tolerates missing grid columns and a single row", {
  results <- data.frame(Comparison = "A", score_sig_network = 0.1)
  expect_equal(drop_constant_grid_columns(results), results)
})

test_that("perc_suffix omits the suffix for 0 and is vectorized", {
  expect_equal(perc_suffix(0), "")
  expect_equal(perc_suffix(0.7), "_0.7")
  expect_equal(perc_suffix(c(0, 0.7)), c("", "_0.7"))
})

test_that("reconstruct/find_missing agree with the actual temp filename when perc_cutoff is 0 (no suffix)", {
  respath <- file.path(tempdir(), "reconstruct_zero_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)
  # Matches what network-score.R's temp_network_file construction actually
  # writes for perc_cutoff = 0: no "_0" suffix.
  writeLines("0.03", file.path(respath, "temp_network_condA.txt"))

  reconstructed <- reconstruct_golden_score_from_temp(respath, conditions = c("condA", "condB"), perc_cutoffs = 0)
  expect_equal(reconstructed$score_sig_network[reconstructed$condition == "condA"], 0.03)

  missing <- find_missing_golden_score_combinations(respath, conditions = c("condA", "condB"), perc_cutoffs = 0)
  expect_equal(missing$condition, "condB")
})
