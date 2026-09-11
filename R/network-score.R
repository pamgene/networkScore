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
#' @param condition_specs A list, one element per grid cell (e.g. built from
#'   `networkGen::build_network_grid()`'s output), each a list with:
#'   `cell_key` (a value unique to this cell -- e.g. every field the grid
#'   varies over, joined -- used to regroup this cell's own observed +
#'   permutation results afterward, so cells that happen to share a
#'   `condition`/`spec_cutoff`/`perc_cutoff` but differ in `b`/`ppi_network`/
#'   `rank_uka_abs` aren't mixed together), `condition` (label),
#'   `spec_cutoff`, `perc_cutoff`, `b`, `rank_uka_abs`, `ppi_network`,
#'   `uka_filt` (observed top-hit kinase data frame), `uka_cell_all` (all
#'   kinase-activity rows for this cell, to build permutations from),
#'   `sens_filt` (observed top-hit sensitivity data frame, or `NULL` for
#'   kinase-only), `vals_extra` (named list merged into this cell's result
#'   row, e.g. `n_kins`), `respath` (folder the observed build is written
#'   into, or `NULL` to not write it -- typically every cell sharing a
#'   combination shares one folder, from `networkGen::prepare_grid_folders()`).
#'   Every spec's `ppi_network` must be the same reference network -- callers
#'   already guarantee this by calling once per output folder, and folders
#'   are keyed by `ppi_network_name` (see `networkGen::prepare_grid_folders()`).
#'   It is read once, from the first spec, and handed to
#'   `networkGen::generate_networks_batch()` as its single shared
#'   `ppi_network` rather than copied into all
#'   `length(condition_specs) * (1 + nPerms)` task bundles -- the latter
#'   serializes a fresh multi-MB copy per task to every `future` worker
#'   (tens of GiB for a real grid; trips `future.globals.maxSize`).
#' @param nPerms Number of permutations per cell.
#' @param generate_fn `networkGen::generate_paired_network` or
#'   `networkGen::generate_kinase_network`.
#' @param uka_top_fn Function `(uka_cell_all, spec_cutoff, rank_uka_abs, perc_cutoff) -> uka_filt`
#'   used to re-filter each permutation's shuffled input -- both kinase-only
#'   and paired paths use `networkGen::uka_top()` here, since by this point
#'   the input has already been cleaned to a common `fscore` column.
#' @param paired If `TRUE`, every task also carries `sens = spec$sens_filt`.
#' @param ... Passed to `generate_fn` for every build (e.g. `n`, `w`, `r`,
#'   `mu`, `seed`) -- uniform across the whole call, not gridded.
#'
#' @return A list, one element per grid cell (same order as
#'   `condition_specs`): `condition`, `vals` (result-row values, including
#'   `obs_network`, `score_sig_network`, `obs_network_inv`,
#'   `score_sig_network_inv`), `metrics_df` (long-format permutation-value
#'   rows), `obs_metrics` (long-format observed-value rows).
#' @keywords internal
score_conditions <- function(condition_specs, nPerms, generate_fn, uka_top_fn, paired = FALSE, ...) {
  extra_args <- list(...)
  # The reference network is constant across the whole call (see @param
  # condition_specs) -- pulled out once here and passed to
  # generate_networks_batch() as its shared `ppi_network`, never placed in
  # a per-task `args`. Putting it in `args` would serialize one copy per
  # task to every worker.
  shared_ppi_network <- if (length(condition_specs) > 0) condition_specs[[1]]$ppi_network else NULL

  build_args <- function(uka, spec, role) {
    write_this <- identical(role, "observed") && !is.null(spec$respath)
    args <- c(extra_args, list(
      uka = uka, condition = spec$condition, spec_cutoff = spec$spec_cutoff,
      b = spec$b, w = spec$w, write = write_this
    ))
    if (write_this) args$res.path <- spec$respath
    if (paired) args$sens <- spec$sens_filt
    args
  }

  tasks <- list()
  for (spec in condition_specs) {
    tasks[[length(tasks) + 1]] <- list(
      args = build_args(spec$uka_filt, spec, role = "observed"),
      meta = list(cell_key = spec$cell_key, condition = spec$condition, spec_cutoff = spec$spec_cutoff, perc_cutoff = spec$perc_cutoff, role = "observed", perm_index = NA_integer_)
    )
    for (i in seq_len(nPerms)) {
      shuffled <- spec$uka_cell_all
      shuffled$uniprotname <- sample(shuffled$uniprotname)
      uka_filt_random <- uka_top_fn(shuffled, spec_cutoff = spec$spec_cutoff, rank_uka_abs = spec$rank_uka_abs, perc_cutoff = spec$perc_cutoff)
      tasks[[length(tasks) + 1]] <- list(
        args = build_args(uka_filt_random, spec, role = "permutation"),
        meta = list(cell_key = spec$cell_key, condition = spec$condition, spec_cutoff = spec$spec_cutoff, perc_cutoff = spec$perc_cutoff, role = "permutation", perm_index = i)
      )
    }
  }

  batch_out <- networkGen::generate_networks_batch(
    tasks, generate_fn = generate_fn, ppi_network = shared_ppi_network
  )

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
    cond_df <- stats_df %>% dplyr::filter(.data$cell_key == spec$cell_key)
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
      dplyr::select("condition", "spec_cutoff", "perc_cutoff", "role", "perm_index", "metric", "value")

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
#' Unified grid interface -- `spec_cutoff`/`perc_cutoff` are gridded (every
#' combination scored, not just paired elementwise), and conditions come
#' from `uka` itself, via [networkGen::build_network_grid()]. No loop over
#' parameter combinations happens here or in [score_conditions()]; the grid
#' is expanded exactly once. The only loop in this function is over the
#' resulting `(spec_cutoff, perc_cutoff)` combinations' *output folders*
#' (one [prepare_score_run_params()] folder each, with its own
#' `results.csv`/checkpointing) -- a consequence of every combination
#' getting its own folder (see `vignette("networkGen")`), not a second grid
#' construction.
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param spec_cutoff,perc_cutoff,b,w,rank_uka_abs Vectors -- every
#'   combination is scored, not just the ones a particular caller happens to
#'   compare. `b`/`w` are the PCSF cost knobs (see
#'   [networkGen::generate_kinase_network()]); both default 2.
#' @param respath Base output directory. Each combination gets its own
#'   parameter-encoded subfolder under it (see [prepare_score_run_params()])
#'   -- `results.csv`, `metrics_permutations.csv`, `metrics_observed.csv`,
#'   the faceted diagnostic histograms, and each condition's observed
#'   network (`nodes_*.csv`/`edges_*.csv`/`wc_df_*.csv`, not permutation
#'   networks) all land there.
#' @param ppi_network A data frame with columns `head`, `tail`, `cost`, or a
#'   fully named list of them (e.g. `list(v12_5 = string_12.5, kins502 =
#'   ppi_networkv12_502_kins)`) to also grid across more than one reference
#'   network.
#' @param nperms_network Number of permutations per condition. Default 30.
#' @param cs `TRUE`/`FALSE` to force per-comparison (csUKA) vs. mean/median
#'   columns; `NULL` (default) auto-detects via [networkGen::detect_csuka()].
#' @param comparison_col Name of the raw UKA column identifying each
#'   comparison, passed through to `networkGen::prep_uka()` and
#'   `networkGen::build_network_grid()`. Different Tercen exports name it
#'   differently -- `"Sgroup_contrast"`, `"Sample"`, ... Default
#'   `"Sgroup_contrast"`.
#' @param max_tasks Refuse to proceed (`stop()`, without building anything)
#'   if the grid, with permutations, expands to more than this many
#'   networks -- a safety guard against an unintentionally huge overnight
#'   run, not a time estimate. Raise explicitly if the size is actually
#'   intended. Default 500.
#' @param ... Passed to every build, uniform across the whole grid (not
#'   gridded) -- e.g. `n`, `r`, `mu`, `seed`.
#'
#' @return A list: `results` (one row per `Comparison` x every grid
#'   combination, across every folder -- `Comparison` is the condition/
#'   contrast label, e.g. `"DrugA_vs_DMSO"`; `spec_cutoff`/`perc_cutoff`/`b`/
#'   `w`/`rank_uka_abs`/`ppi_network` each appear as a column only if that
#'   dimension actually varies across this call's grid, see
#'   [drop_constant_grid_columns()]), `logs`.
#' @export
make_golden_score_kinase <- function(uka, spec_cutoff, perc_cutoff, respath,
                                      ppi_network, b = 2, w = 2, nperms_network = 30,
                                      rank_uka_abs = TRUE, cs = NULL,
                                      comparison_col = "Sgroup_contrast", max_tasks = 500, ...) {
  # Captured immediately, before anything else forces these arguments --
  # forcing a promise before enquo() silently degrades the captured label
  # to a generic value placeholder instead of the caller's actual
  # expression (see networkGen::run_network_grid() for the same fix).
  uka_label <- rlang::as_label(rlang::enquo(uka))
  ppi_network_label <- rlang::as_label(rlang::enquo(ppi_network))

  if (is.null(cs)) cs <- networkGen::detect_csuka(uka)

  grid <- networkGen::build_network_grid(
    uka, clean_fn = function(x) networkGen::prep_uka(x, cs = cs, comparison_col = comparison_col),
    comparison_col = comparison_col,
    spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff, b = b, w = w, rank_uka_abs = rank_uka_abs,
    ppi_network = networkGen::normalize_ppi_network_list(ppi_network, ppi_network_label)
  )
  grid <- purrr::map(grid, function(cell) {
    cell$cell_key <- paste(cell$condition, cell$spec_cutoff, cell$perc_cutoff, cell$b, cell$w, cell$rank_uka_abs, cell$ppi_network_name, sep = "||")
    cell
  })
  cat("Grid has", length(grid), "condition x spec_cutoff x perc_cutoff x b x w x rank_uka_abs x ppi_network cells\n")

  total_tasks <- length(grid) * (1 + nperms_network)
  if (total_tasks > max_tasks) {
    stop(
      "make_golden_score_kinase() would build ", total_tasks, " networks (", length(grid),
      " grid cells x (1 observed + ", nperms_network, " permutations)), over max_tasks = ", max_tasks,
      ". Review the grid before proceeding, or pass a higher max_tasks explicitly if this size ",
      "is actually intended.",
      call. = FALSE
    )
  }

  combos <- networkGen::prepare_grid_folders(
    grid, prepare_fn = prepare_score_run_params, base_labels = list(uka = uka_label),
    respath = respath, uka = uka, nperms_network = nperms_network, cs = cs
  )

  all_results <- list()
  all_logs <- character()

  for (key in names(combos$folders)) {
    folder <- combos$folders[[key]]
    cells <- grid[combos$combo_key == key]

    results <- initialize_or_read_df(folder, "results")
    all_metrics_df <- initialize_or_read_df(folder, "metrics_permutations")
    all_obs_metrics_df <- initialize_or_read_df(folder, "metrics_observed")

    already_done <- if (nrow(results) > 0) results$Comparison else character()
    todo_cells <- Filter(function(cell) !(cell$condition %in% already_done), cells)

    loop_start <- Sys.time()
    temp_files <- c()

    if (length(todo_cells) > 0) {
      condition_specs <- purrr::map(todo_cells, function(cell) {
        list(
          cell_key = cell$cell_key, condition = cell$condition, spec_cutoff = cell$spec_cutoff,
          perc_cutoff = cell$perc_cutoff, b = cell$b, w = cell$w, rank_uka_abs = cell$rank_uka_abs, ppi_network = cell$ppi_network,
          uka_filt = cell$uka_filt, uka_cell_all = cell$uka_cell_all, sens_filt = NULL, respath = folder,
          vals_extra = list(
            Comparison = cell$condition, spec_cutoff = cell$spec_cutoff, perc_cutoff = cell$perc_cutoff,
            b = cell$b, w = cell$w, rank_uka_abs = cell$rank_uka_abs, ppi_network = cell$ppi_network_name,
            n_kins = nrow(cell$uka_filt)
          )
        )
      })

      cond_results <- score_conditions(
        condition_specs, nPerms = nperms_network,
        generate_fn = networkGen::generate_kinase_network,
        uka_top_fn = function(x, spec_cutoff, rank_uka_abs, perc_cutoff) {
          networkGen::uka_top(x, spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
        },
        paired = FALSE, ...
      )

      for (r in cond_results) {
        results <- dplyr::bind_rows(results, as.data.frame(r$vals))
        all_metrics_df <- dplyr::bind_rows(all_metrics_df, r$metrics_df)
        all_obs_metrics_df <- dplyr::bind_rows(all_obs_metrics_df, r$obs_metrics)

        temp_file <- write_network_temp_file(folder, r$condition, r$vals$perc_cutoff, r$vals)
        temp_files <- c(temp_files, temp_file)
      }
    }

    out <- finalize_golden_score_results(results, all_metrics_df, all_obs_metrics_df, temp_files, folder, loop_start, character())
    all_results[[length(all_results) + 1]] <- out$results
    all_logs <- c(all_logs, out$logs)
  }

  list(results = drop_constant_grid_columns(dplyr::bind_rows(all_results)), logs = all_logs)
}

#' Full (paired kinase + sensitivity) golden score analysis
#'
#' Unified grid interface, same shape as [make_golden_score_kinase()]:
#' `spec_cutoff`/`perc_cutoff` are gridded via
#' [networkGen::build_network_grid()], cells are then restricted to cell
#' lines with matching sensitivity data. `sens_filt` depends on `(cell,
#' perc_cutoff)` only (sensitivity has no specificity-cutoff concept), so
#' it's resolved once per cell after the kinase-side grid is built, not
#' gridded independently.
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param sens Raw sensitivity data frame.
#' @param control Control condition name, as it appears in UKA's `contrast` column.
#' @param spec_cutoff,perc_cutoff,b,w,rank_uka_abs Vectors -- every
#'   combination is scored, not just the ones a particular caller happens to
#'   compare. `b`/`w` are the PCSF cost knobs (see
#'   [networkGen::generate_paired_network()]); both default 2.
#' @param respath Base output directory. Each combination gets its own
#'   parameter-encoded subfolder under it (see [prepare_score_run_params()])
#'   -- `results.csv`, `metrics_permutations.csv`, `metrics_observed.csv`,
#'   the faceted diagnostic histograms, and each cell's observed network
#'   (`nodes_*.csv`/`edges_*.csv`/`wc_df_*.csv`, not permutation networks)
#'   all land there.
#' @param uka_fam Kinase-family lookup table (`Kinase_Name`, `Kinase_family`), for overlap scoring.
#' @param ppi_network A data frame with columns `head`, `tail`, `cost`, or a
#'   fully named list of them (e.g. `list(v12_5 = string_12.5, kins502 =
#'   ppi_networkv12_502_kins)`) to also grid across more than one reference
#'   network.
#' @param del_cells Optional character vector of cell lines to exclude.
#' @param zscore If `TRUE`, sensitivity is z-scored across cell lines instead
#'   of fold-change vs. `control`.
#' @param best_drug_per_target Optional data frame restricting sensitivity to
#'   one drug per target.
#' @param score_overlap,score_network If `FALSE`, skip that score entirely
#'   (fields left `NA`). Both default `TRUE`.
#' @param nperms_overlap,nperms_network Number of permutations for each score. Defaults 500, 30.
#' @param balance If `TRUE`, lowers the sensitivity percentile cutoff by 0.2 (see `networkGen::sens_top()`).
#' @param cs `TRUE`/`FALSE` to force per-comparison (csUKA) vs. mean/median
#'   columns; `NULL` (default) auto-detects via [networkGen::detect_csuka()].
#' @param max_tasks Refuse to proceed (`stop()`, without building anything)
#'   if the grid, with permutations, expands to more than this many PCSF
#'   builds -- a safety guard against an unintentionally huge overnight run,
#'   not a time estimate (overlap-score permutations are cheap, no PCSF
#'   calls, and don't count toward this). Raise explicitly if the size is
#'   actually intended. Default 500.
#' @param ... Passed to every network build, uniform across the whole grid
#'   (not gridded) -- e.g. `n`, `r`, `mu`, `seed`.
#'
#' @return A list: `results` (one row per `Comparison` x every grid
#'   combination, across every folder -- `Comparison` is the full,
#'   untruncated "X vs control" contrast label, e.g. `"DrugA vs Control"`;
#'   the cell-line-only join key used internally against sensitivity data
#'   isn't part of the output; `spec_cutoff`/`perc_cutoff`/`b`/`w`/
#'   `rank_uka_abs`/`ppi_network` each appear as a column only if that
#'   dimension actually varies across this call's grid, see
#'   [drop_constant_grid_columns()]), `logs`.
#' @export
make_golden_score_full <- function(uka, sens, control, spec_cutoff, perc_cutoff, respath, uka_fam,
                                    ppi_network, b = 2, w = 2, del_cells = NULL, zscore = FALSE,
                                    best_drug_per_target = NULL, score_overlap = TRUE, score_network = TRUE,
                                    nperms_overlap = 500, nperms_network = 30,
                                    rank_uka_abs = TRUE, balance = FALSE, cs = NULL, max_tasks = 500, ...) {
  # Captured immediately, before anything else forces these arguments --
  # see networkGen::run_network_grid() for why this must happen first.
  uka_label <- rlang::as_label(rlang::enquo(uka))
  sens_label <- rlang::as_label(rlang::enquo(sens))
  ppi_network_label <- rlang::as_label(rlang::enquo(ppi_network))

  if (is.null(cs)) cs <- networkGen::detect_csuka(uka)

  sens_parsed <- prep_sens(sens, control = control, zscore = zscore, best_drug_per_target = best_drug_per_target)

  # spec_cutoff = -Inf disables prep_uka_paired()'s own specificity
  # pre-filter -- uka_top() (inside build_network_grid()) does the real
  # per-cell spec_cutoff filtering instead, since spec_cutoff is now
  # gridded rather than fixed at cleaning time.
  clean_fn <- function(x) prep_uka_paired(x, spec_cutoff = -Inf, control = control, cs = cs)
  grid <- networkGen::build_network_grid(
    uka, clean_fn = clean_fn, comparison_col = "cell_line",
    spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff, b = b, w = w, rank_uka_abs = rank_uka_abs,
    ppi_network = networkGen::normalize_ppi_network_list(ppi_network, ppi_network_label)
  )

  common_cells <- intersect(unique(purrr::map_chr(grid, "condition")), unique(sens_parsed$cell_line))
  if (!is.null(del_cells)) common_cells <- setdiff(common_cells, del_cells)
  grid <- Filter(function(cell) cell$condition %in% common_cells, grid)
  grid <- purrr::map(grid, function(cell) {
    cell$cell_key <- paste(cell$condition, cell$spec_cutoff, cell$perc_cutoff, cell$b, cell$w, cell$rank_uka_abs, cell$ppi_network_name, sep = "||")
    cell
  })
  cat("Grid has", length(grid), "cell x spec_cutoff x perc_cutoff x b x w x rank_uka_abs x ppi_network cells (after matching sensitivity data)\n")

  total_tasks <- length(grid) * (if (score_network) (1 + nperms_network) else 0)
  if (total_tasks > max_tasks) {
    stop(
      "make_golden_score_full() would build ", total_tasks, " networks (", length(grid),
      " grid cells x (1 observed + ", nperms_network, " permutations)), over max_tasks = ", max_tasks,
      ". Review the grid before proceeding, or pass a higher max_tasks explicitly if this size ",
      "is actually intended.",
      call. = FALSE
    )
  }

  combos <- networkGen::prepare_grid_folders(
    grid, prepare_fn = prepare_score_run_params, base_labels = list(uka = uka_label, sens = sens_label),
    respath = respath, uka = uka, sens = sens, nperms_network = nperms_network, cs = cs
  )

  all_results <- list()
  all_logs <- character()

  for (key in names(combos$folders)) {
    folder <- combos$folders[[key]]
    cells <- grid[combos$combo_key == key]

    results <- initialize_or_read_df(folder, "results")
    all_metrics_df <- initialize_or_read_df(folder, "metrics_permutations")
    all_obs_metrics_df <- initialize_or_read_df(folder, "metrics_observed")

    already_done <- if (nrow(results) > 0) results$Comparison else character()
    todo_cells <- Filter(function(cell) !(cell$uka_cell_all$comparison[1] %in% already_done), cells)

    loop_start <- Sys.time()
    temp_files <- c()

    if (length(todo_cells) > 0) {
      cell_sens <- purrr::map(todo_cells, function(cell) {
        sens_parsed %>% dplyr::filter(.data$cell_line == cell$condition) %>% networkGen::sens_top(cell$perc_cutoff, balance = balance)
      })

      network_results <- if (score_network) {
        condition_specs <- purrr::map(seq_along(todo_cells), function(i) {
          cell <- todo_cells[[i]]
          list(
            cell_key = cell$cell_key, condition = cell$condition, spec_cutoff = cell$spec_cutoff,
            perc_cutoff = cell$perc_cutoff, b = cell$b, w = cell$w, rank_uka_abs = cell$rank_uka_abs, ppi_network = cell$ppi_network,
            uka_filt = cell$uka_filt, uka_cell_all = cell$uka_cell_all, sens_filt = cell_sens[[i]],
            respath = folder, vals_extra = list()
          )
        })
        cond_results <- score_conditions(
          condition_specs, nPerms = nperms_network,
          generate_fn = networkGen::generate_paired_network,
          uka_top_fn = function(x, spec_cutoff, rank_uka_abs, perc_cutoff) networkGen::uka_top(x, spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff),
          paired = TRUE, ...
        )
        stats::setNames(cond_results, purrr::map_chr(todo_cells, "cell_key"))
      } else {
        NULL
      }

      for (i in seq_along(todo_cells)) {
        cell <- todo_cells[[i]]
        sens_filt <- cell_sens[[i]]
        vals <- list(
          Comparison = cell$uka_cell_all$comparison[1], spec_cutoff = cell$spec_cutoff, perc_cutoff = cell$perc_cutoff,
          b = cell$b, w = cell$w, rank_uka_abs = cell$rank_uka_abs, ppi_network = cell$ppi_network_name,
          n_targets = nrow(sens_filt), n_kins = nrow(cell$uka_filt), max_sens_value = max(sens_filt$LogFC)
        )

        if (score_overlap) {
          overlap_res <- compute_golden_overlap_scores(
            uka_filt = cell$uka_filt, sens_filt = sens_filt, uka_cell_all = cell$uka_cell_all, uka_fam = uka_fam,
            cell = cell$condition, perc_cutoff = cell$perc_cutoff, nperms_overlap = nperms_overlap, respath = folder,
            spec_cutoff = cell$spec_cutoff, rank_uka_abs = rank_uka_abs
          )
          vals <- c(vals, overlap_res$vals)
          temp_files <- c(temp_files, overlap_res$temp_files)
        } else {
          vals <- c(vals, list(obs_overlap = NA, score_sig_overlap = NA, obs_overlap_fam = NA, score_sig_overlap_fam = NA))
        }

        if (score_network) {
          r <- network_results[[cell$cell_key]]
          vals <- c(vals, r$vals)
          all_metrics_df <- dplyr::bind_rows(all_metrics_df, r$metrics_df)
          all_obs_metrics_df <- dplyr::bind_rows(all_obs_metrics_df, r$obs_metrics)

          temp_network_file <- write_network_temp_file(folder, cell$condition, vals$perc_cutoff, vals)
          temp_files <- c(temp_files, temp_network_file)
        } else {
          vals <- c(vals, list(obs_network = NA, score_sig_network = NA, obs_network_inv = NA, score_sig_network_inv = NA))
        }

        results <- dplyr::bind_rows(results, as.data.frame(vals))
      }
    }

    out <- finalize_golden_score_results(results, all_metrics_df, all_obs_metrics_df, temp_files, folder, loop_start, character())
    all_results[[length(all_results) + 1]] <- out$results
    all_logs <- c(all_logs, out$logs)
  }

  list(results = drop_constant_grid_columns(dplyr::bind_rows(all_results)), logs = all_logs)
}
