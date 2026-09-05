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
    generate_kinase_network = function(uka, condition, spec_cutoff, b, write) {
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

  expect_no_error({
    result <- make_golden_score_kinase(
      uka = uka, spec_cutoff = 0, respath = respath, perc_cutoffs = 0,
      ppi_network = NULL, b = 1, nperms_network = 3
    )
  })

  expect_equal(nrow(result$results), 2)
  expect_true(all(c("score_sig_network", "score_sig_network_inv") %in% colnames(result$results)))
})
