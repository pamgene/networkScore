# Mock generate_fn: builds a "tight" (small, low rel_med_path) network when
# the task's uka carries prize == 1 (used for the observed build in these
# tests), or a uniformly "loose" (larger chain) network otherwise (used for
# every permutation build). This isolates score_conditions()'s flattening/
# regrouping/significance logic from real PCSF and from the specifics of
# percentile-rank filtering (covered separately in test-data-prep.R).
mock_generate_fn <- function(uka, condition, spec_cutoff, b, write, sens = NULL) {
  tight <- isTRUE(uka$prize[1] == 1)
  n <- if (tight) 2 else 12
  edges <- data.frame(from = paste0("N", seq_len(n - 1)), to = paste0("N", seq_len(n - 1) + 1))
  structure(list(
    network = igraph::graph_from_data_frame(edges, directed = FALSE),
    missing_nodes = NULL
  ), class = "networkGen_result")
}

mock_uka_top_fn <- function(x, spec_cutoff, rank_uka_abs, perc_cutoff) {
  data.frame(uniprotname = x$uniprotname, prize = 0) # always "loose", regardless of shuffle
}

make_spec <- function(condition) {
  list(
    condition = condition,
    uka_filt = data.frame(uniprotname = c("K1", "K2"), prize = 1), # "tight" observed
    uka_cell_all = data.frame(uniprotname = paste0("K", 1:6)),
    sens_filt = NULL,
    vals_extra = list(condition = condition, n_kins = 2)
  )
}

test_that("score_conditions flattens all conditions' observed+permutation tasks into results grouped correctly", {
  specs <- list(make_spec("cond_A"), make_spec("cond_B"))
  nPerms <- 5

  results <- score_conditions(
    specs, ppi_network = NULL, spec_cutoff = 0, b = 1, nPerms = nPerms,
    rank_uka_abs = TRUE, perc_cutoff = 0,
    generate_fn = mock_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
  )

  expect_length(results, 2)
  expect_equal(results[[1]]$condition, "cond_A")
  expect_equal(results[[2]]$condition, "cond_B")

  # observed is always the tighter network -> at the floor significance (1/nPerms)
  expect_equal(results[[1]]$vals$score_sig_network, 1 / nPerms)
  expect_equal(results[[2]]$vals$score_sig_network, 1 / nPerms)

  # vals_extra passed through
  expect_equal(results[[1]]$vals$n_kins, 2)

  # exactly nPerms permutation rows and 1 observed row per condition per metric
  expect_equal(nrow(results[[1]]$metrics_df), nPerms * 4) # 4 stats fields
  expect_equal(nrow(results[[1]]$obs_metrics), 1 * 4)
  expect_true(all(results[[1]]$metrics_df$role == "permutation"))
  expect_true(all(results[[1]]$obs_metrics$role == "observed"))
})

test_that("score_conditions returns NA significance when every permutation build fails", {
  failing_generate_fn <- function(uka, condition, spec_cutoff, b, write, sens = NULL) {
    if (isTRUE(uka$prize[1] == 1)) mock_generate_fn(uka, condition, spec_cutoff, b, write, sens) else NULL
  }

  specs <- list(make_spec("cond_A"))
  results <- score_conditions(
    specs, ppi_network = NULL, spec_cutoff = 0, b = 1, nPerms = 3,
    rank_uka_abs = TRUE, perc_cutoff = 0,
    generate_fn = failing_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
  )

  expect_true(is.na(results[[1]]$vals$score_sig_network))
})

test_that("score_conditions includes sens in build args when paired = TRUE", {
  captured_args <- list()
  capturing_generate_fn <- function(uka, sens, condition, spec_cutoff, b, write) {
    captured_args[[length(captured_args) + 1]] <<- sens
    mock_generate_fn(uka, condition, spec_cutoff, b, write)
  }

  spec <- make_spec("cond_A")
  spec$sens_filt <- data.frame(name = "S1", prize = 0.5)

  score_conditions(
    list(spec), ppi_network = NULL, spec_cutoff = 0, b = 1, nPerms = 1,
    rank_uka_abs = TRUE, perc_cutoff = 0,
    generate_fn = capturing_generate_fn, uka_top_fn = mock_uka_top_fn, paired = TRUE
  )

  expect_true(length(captured_args) >= 1)
  expect_equal(captured_args[[1]], spec$sens_filt)
})
