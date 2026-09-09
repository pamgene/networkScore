# Regression test for a real bug found via end-to-end testing with real PCSF:
# uka_top_fn (used to re-filter each permutation's shuffled input) was being
# called on already-cleaned data (clean_uka_to_kinograte_kinase()'s output,
# columns Sample/uniprotname/LogFC/fscore), but the old uka_top_kinase()
# filtered on the raw "Mean Specificity Score"/"Specificity Score" column
# name -- which no longer exists post-cleaning. Fixed by using
# networkGen::uka_top() (which filters on `fscore`) throughout instead.
# This test exercises the full clean -> filter -> (mocked) generate chain
# without needing real PCSF, so this class of column-name mismatch gets
# caught by the regular test suite going forward.
test_that("make_golden_score_kinase's permutation filtering uses the same column contract as the cleaned data it's given", {
  local_mocked_bindings(
    generate_kinase_network = function(uka, condition, spec_cutoff, b, ppi_network = NULL, write, res.path = NULL) {
      structure(list(
        network = igraph::graph_from_data_frame(data.frame(from = uka$name[1], to = uka$name[min(2, nrow(uka))]), directed = FALSE),
        missing_nodes = NULL
      ), class = "networkGen_result")
    },
    .package = "networkGen"
  )

  uka <- data.frame(
    `x.Sample` = rep(c("cond_A", "cond_B"), each = 3),
    `x.Kinase Name` = rep(c("K1", "K2", "K3"), 2),
    `x.Median Kinase Statistic` = c(2.0, -1.5, 0.3, 1.8, -1.2, 0.1),
    `x.Mean Specificity Score` = 2.0,
    check.names = FALSE
  )

  respath <- file.path(tempdir(), "integration_kinase_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  # Named variable, not a literal -- as_label() on a literal data.frame() call
  # would embed its quote/paren characters into the output folder name.
  fake_ppi <- data.frame(head = "A", tail = "B", cost = 0.1)

  expect_no_error({
    result <- make_golden_score_kinase(
      uka = uka, spec_cutoff = 0, respath = respath, perc_cutoff = 0,
      ppi_network = fake_ppi, b = 1, nperms_network = 3
    )
  })

  expect_equal(nrow(result$results), 2)
  expect_true(all(c("score_sig_network", "score_sig_network_inv") %in% colnames(result$results)))
})

# Regression test for Candidate D (score_conditions() caller duplication,
# see DevOpti/networkGen/docs/adr and the codebase-architecture review):
# make_golden_score_full() used to build a one-cell condition_specs list
# inside its per-cell loop and call score_conditions() once per cell,
# defeating full flattening (ADR 0003) for the paired path even though
# make_golden_score_kinase() already did this correctly. Fixed by building
# every todo cell's inputs up front and submitting one score_conditions()
# call spanning all of them, per perc_cutoff -- this test locks that in by
# counting generate_networks_batch() calls and each call's task-list size
# directly, rather than only checking the final results (which would look
# identical either way).
test_that("make_golden_score_full batches all cells' network-score builds into a single generate_networks_batch call per perc_cutoff (full flattening)", {
  cleaned_uka <- data.frame(
    cell_line = rep(c("cellA", "cellB"), each = 3),
    comparison = rep(c("cellA vs DMSO", "cellB vs DMSO"), each = 3),
    uniprotname = rep(c("K1", "K2", "K3"), 2),
    LogFC = c(2.0, -1.5, 0.3, 1.8, -1.2, 0.1),
    fscore = 2.0
  )
  cleaned_sens <- data.frame(
    cell_line = rep(c("cellA", "cellB"), each = 2),
    uniprotname = rep(c("T1", "T2"), 2),
    LogFC = c(1.0, 0.5, 1.0, 0.5)
  )

  local_mocked_bindings(
    clean_uka_to_kinograte_full = function(uka, spec_cutoff, control, cs = FALSE) cleaned_uka,
    clean_sens_to_kinograte = function(sens, control, zscore = FALSE, del_cell = NULL, best_drug_per_target = NULL) cleaned_sens,
    .package = "networkScore"
  )

  batch_calls <- list()
  local_mocked_bindings(
    generate_networks_batch = function(tasks, generate_fn, ppi_network = NULL, extra_args = list(), progress = FALSE) {
      batch_calls[[length(batch_calls) + 1]] <<- length(tasks)
      lapply(tasks, function(task) {
        tight <- identical(task$meta$role, "observed")
        n <- if (tight) 2 else 12
        edges <- data.frame(from = paste0("N", seq_len(n - 1)), to = paste0("N", seq_len(n - 1) + 1))
        list(
          result = structure(list(
            network = igraph::graph_from_data_frame(edges, directed = FALSE),
            missing_nodes = NULL
          ), class = "networkGen_result"),
          meta = task$meta
        )
      })
    },
    .package = "networkGen"
  )

  respath <- file.path(tempdir(), "integration_full_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  # Real callers always pass variables here (data frames), never string
  # literals -- as_label() on a literal includes its quote characters,
  # which are illegal in a file path. Mirror real usage in the fixture.
  raw_uka_unused <- "raw_uka_unused"
  raw_sens_unused <- "raw_sens_unused"
  fake_ppi <- data.frame(head = "A", tail = "B", cost = 0.1)

  result <- make_golden_score_full(
    uka = raw_uka_unused, sens = raw_sens_unused, control = "DMSO",
    spec_cutoff = 0, respath = respath, uka_fam = NULL, perc_cutoff = 0,
    ppi_network = fake_ppi, b = 1, nperms_network = 3,
    score_overlap = FALSE, score_network = TRUE
  )

  # 2 cells x (1 observed + 3 permutations) = 8 tasks, in ONE batch call --
  # not two separate 4-task calls (the pre-fix behavior).
  expect_equal(batch_calls, list(8))

  expect_equal(nrow(result$results), 2)
  expect_setequal(result$results$Comparison, c("cellA vs DMSO", "cellB vs DMSO"))
  expect_true(all(c("score_sig_network", "score_sig_network_inv") %in% colnames(result$results)))
})

test_that("make_golden_score_kinase gives each (spec_cutoff, perc_cutoff) combination its own output folder", {
  local_mocked_bindings(
    generate_kinase_network = function(uka, condition, spec_cutoff, b, ppi_network = NULL, write, res.path = NULL) {
      structure(list(
        network = igraph::graph_from_data_frame(data.frame(from = "K1", to = "K2"), directed = FALSE),
        missing_nodes = NULL
      ), class = "networkGen_result")
    },
    .package = "networkGen"
  )

  uka <- data.frame(
    `x.Sample` = rep("cond_A", 3),
    `x.Kinase Name` = c("K1", "K2", "K3"),
    `x.Median Kinase Statistic` = c(2.0, -1.5, 0.3),
    `x.Mean Specificity Score` = 2.0,
    check.names = FALSE
  )

  respath <- file.path(tempdir(), "integration_kinase_multicombo_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)
  fake_ppi <- data.frame(head = "A", tail = "B", cost = 0.1)

  result <- make_golden_score_kinase(
    uka = uka, spec_cutoff = c(0, 0.5), respath = respath, perc_cutoff = 0,
    ppi_network = fake_ppi, b = 1, nperms_network = 2
  )

  # 2 spec_cutoff values -> 2 separate result-writing folders
  result_files <- list.files(respath, pattern = "^results\\.csv$", recursive = TRUE, full.names = TRUE)
  expect_length(result_files, 2)
  expect_equal(nrow(result$results), 2) # 1 condition x 2 spec_cutoff
})

test_that("make_golden_score_kinase refuses to proceed when the grid (with permutations) exceeds max_tasks", {
  uka <- data.frame(
    `x.Sample` = rep(paste0("cond", 1:5), each = 2),
    `x.Kinase Name` = rep(c("K1", "K2"), 5),
    `x.Median Kinase Statistic` = 1,
    `x.Mean Specificity Score` = 2.0,
    check.names = FALSE
  )

  respath <- file.path(tempdir(), "integration_kinase_maxtasks_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  expect_error(
    make_golden_score_kinase(
      uka = uka, spec_cutoff = 0, respath = respath, perc_cutoff = 0,
      ppi_network = data.frame(head = "A", tail = "B", cost = 0.1), b = 1, nperms_network = 50, max_tasks = 10
    ),
    "max_tasks"
  )
})
