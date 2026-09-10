# Mock generate_fn: builds a "tight" (small, low rel_med_path) network when
# the task's uka carries prize == 1 (used for the observed build in these
# tests), or a uniformly "loose" (larger chain) network otherwise (used for
# every permutation build). This isolates score_conditions()'s flattening/
# regrouping/significance logic from real PCSF and from the specifics of
# percentile-rank filtering (covered separately in test-data-prep.R).
mock_generate_fn <- function(uka, condition, spec_cutoff, b, w = NULL, ppi_network = NULL, write, sens = NULL, ...) {
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

make_spec <- function(condition, spec_cutoff = 0, perc_cutoff = 0, respath = NULL, b = 1, w = 2, rank_uka_abs = TRUE, ppi_network = NULL) {
  list(
    cell_key = paste(condition, spec_cutoff, perc_cutoff, b, w, rank_uka_abs, sep = "||"),
    condition = condition, spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff, b = b, w = w,
    rank_uka_abs = rank_uka_abs, ppi_network = ppi_network, respath = respath,
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
    specs, nPerms = nPerms,
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

test_that("score_conditions keeps cells with the same condition but different spec_cutoff/perc_cutoff distinct", {
  specs <- list(
    make_spec("cond_A", spec_cutoff = 0, perc_cutoff = 0),
    make_spec("cond_A", spec_cutoff = 0.5, perc_cutoff = 0)
  )

  results <- score_conditions(
    specs, nPerms = 3,
    generate_fn = mock_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
  )

  expect_length(results, 2)
  expect_equal(purrr::map_dbl(results, ~ .x$vals$score_sig_network), c(1 / 3, 1 / 3))
})

test_that("score_conditions returns NA significance when every permutation build fails", {
  failing_generate_fn <- function(uka, condition, spec_cutoff, b, w = NULL, ppi_network = NULL, write, sens = NULL, ...) {
    if (isTRUE(uka$prize[1] == 1)) mock_generate_fn(uka, condition, spec_cutoff, b, w, ppi_network, write, sens) else NULL
  }

  specs <- list(make_spec("cond_A"))
  results <- score_conditions(
    specs, nPerms = 3,
    generate_fn = failing_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
  )

  expect_true(is.na(results[[1]]$vals$score_sig_network))
})

test_that("score_conditions writes the observed build (write=TRUE, res.path=spec$respath) but never a permutation build", {
  captured_calls <- list()
  capturing_generate_fn <- function(uka, condition, spec_cutoff, b, w = NULL, ppi_network = NULL, write, res.path = NULL, sens = NULL, ...) {
    captured_calls[[length(captured_calls) + 1]] <<- list(tight = isTRUE(uka$prize[1] == 1), write = write, res.path = res.path)
    mock_generate_fn(uka, condition, spec_cutoff, b, w, ppi_network, write, sens)
  }

  score_conditions(
    list(make_spec("cond_A", respath = "some/output/path")), nPerms = 3,
    generate_fn = capturing_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
  )

  observed_calls <- Filter(function(c) c$tight, captured_calls)
  perm_calls <- Filter(function(c) !c$tight, captured_calls)

  expect_length(observed_calls, 1)
  expect_true(observed_calls[[1]]$write)
  expect_equal(observed_calls[[1]]$res.path, "some/output/path")

  expect_length(perm_calls, 3)
  expect_true(all(!vapply(perm_calls, `[[`, logical(1), "write")))
  expect_true(all(vapply(perm_calls, function(c) is.null(c$res.path), logical(1))))
})

test_that("score_conditions never writes any build when respath is not given on the spec (default, matches previous behavior)", {
  captured_writes <- c()
  capturing_generate_fn <- function(uka, condition, spec_cutoff, b, w = NULL, ppi_network = NULL, write, res.path = NULL, sens = NULL, ...) {
    captured_writes <<- c(captured_writes, write)
    mock_generate_fn(uka, condition, spec_cutoff, b, w, ppi_network, write, sens)
  }

  score_conditions(
    list(make_spec("cond_A")), nPerms = 3,
    generate_fn = capturing_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
  )

  expect_true(all(!captured_writes))
})

test_that("score_conditions hands one shared ppi_network to the batch, not a per-task copy", {
  shared <- data.frame(head = "A", tail = "B", cost = 0.1)
  specs <- list(
    make_spec("cond_A", ppi_network = shared),
    make_spec("cond_B", ppi_network = shared)
  )

  captured <- NULL
  testthat::with_mocked_bindings(
    generate_networks_batch = function(tasks, generate_fn, ppi_network = NULL, extra_args = list(), progress = FALSE) {
      captured <<- list(
        n_tasks = length(tasks),
        any_task_carries_ppi = any(vapply(tasks, function(t) "ppi_network" %in% names(t$args), logical(1))),
        shared_ppi = ppi_network
      )
      lapply(tasks, function(t) list(result = NULL, meta = t$meta))
    },
    .package = "networkGen",
    {
      score_conditions(
        specs, nPerms = 2,
        generate_fn = mock_generate_fn, uka_top_fn = mock_uka_top_fn, paired = FALSE
      )
    }
  )

  # 2 conditions x (1 observed + 2 permutations) = 6 tasks, none of which
  # embeds the multi-MB reference network in its args; the single shared
  # network is passed once as generate_networks_batch()'s `ppi_network`.
  expect_equal(captured$n_tasks, 6)
  expect_false(captured$any_task_carries_ppi)
  expect_identical(captured$shared_ppi, shared)
})

test_that("score_conditions includes sens in build args when paired = TRUE", {
  captured_args <- list()
  capturing_generate_fn <- function(uka, sens, condition, spec_cutoff, b, w = NULL, ppi_network = NULL, write, ...) {
    captured_args[[length(captured_args) + 1]] <<- sens
    mock_generate_fn(uka, condition, spec_cutoff, b, w, ppi_network, write)
  }

  spec <- make_spec("cond_A")
  spec$sens_filt <- data.frame(name = "S1", prize = 0.5)

  score_conditions(
    list(spec), nPerms = 1,
    generate_fn = capturing_generate_fn, uka_top_fn = mock_uka_top_fn, paired = TRUE
  )

  expect_true(length(captured_args) >= 1)
  expect_equal(captured_args[[1]], spec$sens_filt)
})
