#' Score every condition's network significance in one flattened batch
#'
#' The core full-flattening engine (see `networkGen`'s
#' `docs/adr/0003-full-flattening-parallelization.md` -- suite-wide ADRs are
#' filed in `networkGen`'s repo even when, as here, the decision is about
#' `networkScore`'s own design): builds one flat list of every network build needed
#' -- one observed + `nPerms` permutation builds, for every condition in
#' `condition_specs` -- submits it as a single
#' `networkGen::generate_networks_batch()` call, then regroups the results by
#' condition to compute each condition's significance score. Not exported;
#' called by [make_golden_score_full()] / [make_golden_score_kinase()].
#'
#' @param condition_specs A list, one element per condition, each a list
#'   with: `condition` (label), `uka_filt` (observed top-hit kinase data
#'   frame), `uka_cell_all` (all kinase-activity rows for this condition, to
#'   build permutations from), `sens_filt` (observed top-hit sensitivity
#'   data frame, or `NULL` for kinase-only), `vals_extra` (named list merged
#'   into this condition's result row, e.g. `n_kins`).
#' @param ppi_network Data frame with columns `head`, `tail`, `cost`.
#' @param spec_cutoff,b Passed to the generate function for every build.
#' @param nPerms Number of permutations per condition.
#' @param rank_uka_abs,perc_cutoff Passed to `uka_top_fn` when building each
#'   permutation's shuffled top-hit input.
#' @param generate_fn `networkGen::generate_paired_network` or
#'   `networkGen::generate_kinase_network`.
#' @param uka_top_fn Function `(uka_cell_all, spec_cutoff, rank_uka_abs, perc_cutoff) -> uka_filt`
#'   used to re-filter each permutation's shuffled input -- both kinase-only
#'   and paired paths use `networkGen::uka_top()` here, since by this point
#'   the input has already been cleaned to a common `fscore` column.
#' @param paired If `TRUE`, every task also carries `sens = spec$sens_filt`.
#'
#' @return A list, one element per condition (same order as `condition_specs`):
#'   `condition`, `vals` (result-row values, including `obs_network`,
#'   `score_sig_network`, `obs_network_inv`, `score_sig_network_inv`),
#'   `metrics_df` (long-format permutation-value rows), `obs_metrics`
#'   (long-format observed-value rows).
#' @param respath If given, the observed build for each condition is written
#'   (`write = TRUE`, `res.path = respath`) via [networkGen::kinograte_pg_pcsf()]'s
#'   usual `nodes_*.csv`/`edges_*.csv`/`wc_df_*.csv` output -- permutation
#'   builds are never written regardless (there are typically dozens per
#'   condition and they're not independently meaningful). `NULL` (default):
#'   no build is written, matching the previous behavior.
#' @keywords internal
score_conditions <- function(condition_specs, ppi_network, spec_cutoff, b, nPerms,
                              rank_uka_abs, perc_cutoff, generate_fn, uka_top_fn, paired = FALSE,
                              respath = NULL) {
  build_args <- function(uka, condition, sens, role) {
    write_this <- identical(role, "observed") && !is.null(respath)
    args <- list(uka = uka, condition = condition, spec_cutoff = spec_cutoff, b = b, write = write_this)
    if (write_this) args$res.path <- respath
    if (paired) args$sens <- sens
    args
  }

  tasks <- list()
  for (spec in condition_specs) {
    tasks[[length(tasks) + 1]] <- list(
      args = build_args(spec$uka_filt, spec$condition, spec$sens_filt, role = "observed"),
      meta = list(condition = spec$condition, role = "observed", perm_index = NA_integer_)
    )
    for (i in seq_len(nPerms)) {
      shuffled <- spec$uka_cell_all
      shuffled$uniprotname <- sample(shuffled$uniprotname)
      uka_filt_random <- uka_top_fn(shuffled, spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      tasks[[length(tasks) + 1]] <- list(
        args = build_args(uka_filt_random, spec$condition, spec$sens_filt, role = "permutation"),
        meta = list(condition = spec$condition, role = "permutation", perm_index = i)
      )
    }
  }

  batch_out <- networkGen::generate_networks_batch(tasks, generate_fn = generate_fn, ppi_network = ppi_network)

  stats_rows <- purrr::map(batch_out, function(item) {
    s <- if (is.null(item$result)) {
      list(rel_med_path = NA_real_, rel_med_path_inv = NA_real_, density = NA_real_, clustering = NA_real_)
    } else {
      compute_network_stats(item$result)
    }
    c(item$meta, s)
  })
  stats_df <- dplyr::bind_rows(stats_rows)

  purrr::map(condition_specs, function(spec) {
    cond_df <- stats_df %>% dplyr::filter(.data$condition == spec$condition)
    obs_row <- cond_df %>% dplyr::filter(.data$role == "observed")
    perm_rows <- cond_df %>% dplyr::filter(.data$role == "permutation")

    score_obs <- obs_row$rel_med_path[1]
    score_obs_inv <- obs_row$rel_med_path_inv[1]
    nPerms_notna <- sum(!is.na(perm_rows$rel_med_path))

    if (is.na(score_obs) || nPerms_notna == 0) {
      vals <- c(spec$vals_extra, list(
        obs_network = NA_real_, score_sig_network = NA_real_,
        obs_network_inv = NA_real_, score_sig_network_inv = NA_real_
      ))
    } else {
      sig <- max(sum(perm_rows$rel_med_path <= score_obs, na.rm = TRUE) / nPerms_notna, 1 / nPerms_notna)
      sig_inv <- max(sum(perm_rows$rel_med_path_inv >= score_obs_inv, na.rm = TRUE) / nPerms_notna, 1 / nPerms_notna)
      vals <- c(spec$vals_extra, list(
        obs_network = round(score_obs, 3), score_sig_network = sig,
        obs_network_inv = round(score_obs_inv, 3), score_sig_network_inv = sig_inv
      ))
    }

    metrics_long <- cond_df %>%
      tidyr::pivot_longer(
        cols = c("rel_med_path", "rel_med_path_inv", "density", "clustering"),
        names_to = "metric", values_to = "value"
      ) %>%
      dplyr::mutate(perc_cutoff = perc_cutoff) %>%
      dplyr::select("condition", "perc_cutoff", "role", "perm_index", "metric", "value")

    list(
      condition = spec$condition, vals = vals,
      metrics_df = metrics_long %>% dplyr::filter(.data$role == "permutation"),
      obs_metrics = metrics_long %>% dplyr::filter(.data$role == "observed")
    )
  })
}

#' Write a condition's network-score temp progress file
#'
#' Shared by [make_golden_score_kinase()] / [make_golden_score_full()],
#' called once per condition right after [score_conditions()] returns.
#' Mirrors the `temp_overlap_*` files [compute_golden_overlap_scores()]
#' writes -- deleted by [finalize_golden_score_results()] once the whole run
#' completes; lets a crashed/interrupted run be resumed.
#'
#' @param respath Folder to write into.
#' @param condition Condition/cell-line label.
#' @param perc_cutoff Percentile cutoff, for the file name (via [perc_suffix()]).
#' @param vals A condition's result-row values, from `score_conditions()`'s
#'   per-condition `vals` (must include `score_sig_network`/`score_sig_network_inv`).
#'
#' @return The path written, invisibly.
#' @keywords internal
write_network_temp_file <- function(respath, condition, perc_cutoff, vals) {
  temp_file <- file.path(respath, paste0("temp_network_", condition, perc_suffix(perc_cutoff), ".txt"))
  writeLines(c(as.character(vals$score_sig_network), as.character(vals$score_sig_network_inv)), temp_file)
  invisible(temp_file)
}

#' Unified golden-score entry point
#'
#' Dispatches to [make_golden_score_full()] (paired kinase+sensitivity) when
#' `sens` is given, or [make_golden_score_kinase()] (kinase-only) otherwise.
#'
#' @param uka UKA (kinase-activity) data.
#' @param sens Optional sensitivity data; if `NULL`, runs kinase-only mode.
#' @param ... Passed through to the dispatched function.
#'
#' @return See [make_golden_score_full()] / [make_golden_score_kinase()].
#' @export
make_golden_score <- function(uka, sens = NULL, ...) {
  if (is.null(sens)) {
    cat("Running kinase-only golden score analysis...\n")
    make_golden_score_kinase(uka = uka, ...)
  } else {
    cat("Running full golden score analysis (UKA + sensitivity)...\n")
    make_golden_score_full(uka = uka, sens = sens, ...)
  }
}

#' Kinase-only golden score analysis
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param spec_cutoff Specificity-score cutoff.
#' @param respath Base output directory. The actual results go into a
#'   parameter-encoded subfolder under it (see [prepare_score_run_params()]) --
#'   `results.csv`, `metrics_permutations.csv`, `metrics_observed.csv`, the
#'   faceted diagnostic histograms, and each condition's observed network
#'   (`nodes_*.csv`/`edges_*.csv`/`wc_df_*.csv`, not permutation networks)
#'   all land there.
#' @param perc_cutoffs Vector of percentile cutoffs to score at.
#' @param ppi_network Data frame with columns `head`, `tail`, `cost`.
#' @param b PCSF terminal-prize weight.
#' @param nperms_network Number of permutations per condition. Default 50.
#' @param rank_uka_abs If `TRUE` (default), rank kinase hits by `abs(LogFC)`.
#' @param cs If `TRUE`, use the per-comparison specificity column.
#'
#' @return A list: `results` (one row per condition x perc_cutoff), `logs`.
#' @export
make_golden_score_kinase <- function(uka, spec_cutoff, respath, perc_cutoffs,
                                      ppi_network, b, nperms_network = 50,
                                      rank_uka_abs = TRUE, cs = FALSE) {
  run_params <- prepare_score_run_params(
    respath = respath, uka = uka, ppi_network = ppi_network, spec_cutoff = spec_cutoff,
    b = b, nperms_network = nperms_network, rank_uka_abs = rank_uka_abs, cs = cs
  )
  respath <- run_params$respath

  uka_parsed <- clean_uka_to_kinograte_kinase(uka, cs = cs)
  conditions <- unique(uka_parsed$Sample)
  cat("Dataset has", length(conditions), "conditions for processing\n")

  results <- initialize_or_read_df(respath, "results")
  all_metrics_df <- initialize_or_read_df(respath, "metrics_permutations")
  all_obs_metrics_df <- initialize_or_read_df(respath, "metrics_observed")

  loop_start <- Sys.time()
  logs <- character()
  temp_files <- c()

  for (perc_cutoff in perc_cutoffs) {
    already_done <- if (nrow(results) > 0) {
      conditions[conditions %in% results$condition[results$perc_cutoff == perc_cutoff]]
    } else {
      character()
    }
    todo <- setdiff(conditions, already_done)
    if (length(todo) == 0) next

    condition_specs <- purrr::map(todo, function(condition) {
      uka_filt <- uka_parsed %>%
        dplyr::filter(.data$Sample == condition) %>%
        networkGen::uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      uka_cell_all <- uka_parsed %>% dplyr::filter(.data$Sample == condition)
      list(
        condition = condition, uka_filt = uka_filt, uka_cell_all = uka_cell_all, sens_filt = NULL,
        vals_extra = list(condition = condition, spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff, n_kins = nrow(uka_filt))
      )
    })

    cond_results <- score_conditions(
      condition_specs, ppi_network = ppi_network, spec_cutoff = spec_cutoff, b = b,
      nPerms = nperms_network, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff,
      generate_fn = networkGen::generate_kinase_network,
      uka_top_fn = function(x, spec_cutoff, rank_uka_abs, perc_cutoff) {
        networkGen::uka_top(x, spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      },
      paired = FALSE, respath = respath
    )

    for (r in cond_results) {
      results <- dplyr::bind_rows(results, as.data.frame(r$vals))
      all_metrics_df <- dplyr::bind_rows(all_metrics_df, r$metrics_df)
      all_obs_metrics_df <- dplyr::bind_rows(all_obs_metrics_df, r$obs_metrics)

      temp_file <- write_network_temp_file(respath, r$condition, perc_cutoff, r$vals)
      temp_files <- c(temp_files, temp_file)
    }
  }

  finalize_golden_score_results(results, all_metrics_df, all_obs_metrics_df, temp_files, respath, loop_start, logs)
}

#' Full (paired kinase + sensitivity) golden score analysis
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param sens Raw sensitivity data frame.
#' @param control Control condition name, as it appears in UKA's `contrast` column.
#' @param spec_cutoff Specificity-score cutoff.
#' @param respath Base output directory. The actual results go into a
#'   parameter-encoded subfolder under it (see [prepare_score_run_params()]) --
#'   `results.csv`, `metrics_permutations.csv`, `metrics_observed.csv`, the
#'   faceted diagnostic histograms, and each cell's observed network
#'   (`nodes_*.csv`/`edges_*.csv`/`wc_df_*.csv`, not permutation networks)
#'   all land there.
#' @param uka_fam Kinase-family lookup table (`Kinase_Name`, `Kinase_family`), for overlap scoring.
#' @param perc_cutoffs Vector of percentile cutoffs to score at.
#' @param ppi_network Data frame with columns `head`, `tail`, `cost`.
#' @param b PCSF terminal-prize weight.
#' @param del_cells Optional character vector of cell lines to exclude.
#' @param zscore If `TRUE`, sensitivity is z-scored across cell lines instead
#'   of fold-change vs. `control`.
#' @param best_drug_per_target Optional data frame restricting sensitivity to
#'   one drug per target.
#' @param score_overlap,score_network If `FALSE`, skip that score entirely
#'   (fields left `NA`). Both default `TRUE`.
#' @param nperms_overlap,nperms_network Number of permutations for each score. Defaults 500, 50.
#' @param rank_uka_abs If `TRUE` (default), rank kinase hits by `abs(LogFC)`.
#' @param balance If `TRUE`, lowers the sensitivity percentile cutoff by 0.2 (see `networkGen::sens_top()`).
#' @param cs If `TRUE`, use the per-comparison specificity column.
#'
#' @return A list: `results` (one row per cell x perc_cutoff), `logs`.
#' @export
make_golden_score_full <- function(uka, sens, control, spec_cutoff, respath, uka_fam, perc_cutoffs,
                                    ppi_network, b, del_cells = NULL, zscore = FALSE,
                                    best_drug_per_target = NULL, score_overlap = TRUE, score_network = TRUE,
                                    nperms_overlap = 500, nperms_network = 50,
                                    rank_uka_abs = TRUE, balance = FALSE, cs = FALSE) {
  run_params <- prepare_score_run_params(
    respath = respath, uka = uka, sens = sens, ppi_network = ppi_network, spec_cutoff = spec_cutoff,
    b = b, nperms_network = nperms_network, rank_uka_abs = rank_uka_abs, cs = cs
  )
  respath <- run_params$respath

  sens_parsed <- clean_sens_to_kinograte(sens, control = control, zscore = zscore, best_drug_per_target = best_drug_per_target)
  uka_parsed <- clean_uka_to_kinograte_full(uka, spec_cutoff = spec_cutoff, control = control, cs = cs)

  common_cells <- intersect(unique(uka_parsed$cell_line), unique(sens_parsed$cell_line))
  if (!is.null(del_cells)) common_cells <- setdiff(common_cells, del_cells)

  results <- initialize_or_read_df(respath, "results")
  all_metrics_df <- initialize_or_read_df(respath, "metrics_permutations")
  all_obs_metrics_df <- initialize_or_read_df(respath, "metrics_observed")

  loop_start <- Sys.time()
  logs <- character()
  temp_files <- c()

  for (perc_cutoff in perc_cutoffs) {
    already_done <- if (nrow(results) > 0) {
      common_cells[common_cells %in% results$cell[results$perc_cutoff == perc_cutoff]]
    } else {
      character()
    }
    todo <- setdiff(common_cells, already_done)
    if (length(todo) == 0) next

    # Prepare every todo cell's filtered inputs up front, same as
    # make_golden_score_kinase(), so the network-score path below can submit
    # one score_conditions() call spanning all cells instead of one per cell
    # -- see docs/adr/0003-full-flattening-parallelization.md in networkGen.
    cell_inputs <- purrr::map(todo, function(cell) {
      uka_filt <- uka_parsed %>% dplyr::filter(.data$cell_line == cell) %>% networkGen::uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      sens_filt <- sens_parsed %>% dplyr::filter(.data$cell_line == cell) %>% networkGen::sens_top(perc_cutoff, balance = balance)
      uka_cell_all <- uka_parsed %>% dplyr::filter(.data$cell_line == cell)
      list(cell = cell, uka_filt = uka_filt, sens_filt = sens_filt, uka_cell_all = uka_cell_all)
    })
    names(cell_inputs) <- todo

    network_results <- if (score_network) {
      condition_specs <- purrr::map(cell_inputs, function(ci) {
        list(condition = ci$cell, uka_filt = ci$uka_filt, uka_cell_all = ci$uka_cell_all, sens_filt = ci$sens_filt, vals_extra = list())
      })
      cond_results <- score_conditions(
        condition_specs, ppi_network = ppi_network, spec_cutoff = spec_cutoff, b = b, nPerms = nperms_network,
        rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff,
        generate_fn = networkGen::generate_paired_network,
        uka_top_fn = function(x, spec_cutoff, rank_uka_abs, perc_cutoff) networkGen::uka_top(x, spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff),
        paired = TRUE, respath = respath
      )
      stats::setNames(cond_results, todo)
    } else {
      NULL
    }

    for (cell in todo) {
      ci <- cell_inputs[[cell]]
      vals <- list(cell = cell, perc_cutoff = perc_cutoff, n_targets = nrow(ci$sens_filt), n_kins = nrow(ci$uka_filt), max_sens_value = max(ci$sens_filt$LogFC))

      if (score_overlap) {
        overlap_res <- compute_golden_overlap_scores(
          uka_filt = ci$uka_filt, sens_filt = ci$sens_filt, uka_cell_all = ci$uka_cell_all, uka_fam = uka_fam,
          cell = cell, perc_cutoff = perc_cutoff, nperms_overlap = nperms_overlap, respath = respath,
          spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs
        )
        vals <- c(vals, overlap_res$vals)
        temp_files <- c(temp_files, overlap_res$temp_files)
      } else {
        vals <- c(vals, list(obs_overlap = NA, score_sig_overlap = NA, obs_overlap_fam = NA, score_sig_overlap_fam = NA))
      }

      if (score_network) {
        r <- network_results[[cell]]
        vals <- c(vals, r$vals)
        all_metrics_df <- dplyr::bind_rows(all_metrics_df, r$metrics_df)
        all_obs_metrics_df <- dplyr::bind_rows(all_obs_metrics_df, r$obs_metrics)

        temp_network_file <- write_network_temp_file(respath, cell, perc_cutoff, vals)
        temp_files <- c(temp_files, temp_network_file)
      } else {
        vals <- c(vals, list(obs_network = NA, score_sig_network = NA, obs_network_inv = NA, score_sig_network_inv = NA))
      }

      results <- dplyr::bind_rows(results, as.data.frame(vals))
    }
  }

  finalize_golden_score_results(results, all_metrics_df, all_obs_metrics_df, temp_files, respath, loop_start, logs)
}
