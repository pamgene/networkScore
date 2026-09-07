#' Read an existing results data frame, or start a fresh empty one
#'
#' @param respath Results folder.
#' @param dfname Base file name (without `.csv`), e.g. `"results"`.
#'
#' @return A data frame -- the existing `<dfname>.csv` if present (so a
#'   resumed run adds to it), or an empty `data.frame()` otherwise.
#' @export
initialize_or_read_df <- function(respath, dfname) {
  df_filename <- paste0(respath, "/", dfname, ".csv")
  if (file.exists(df_filename)) {
    print(paste0(dfname, " df already exists. New data will be added to the existing one."))
    readr::read_csv(df_filename, show_col_types = FALSE)
  } else {
    data.frame()
  }
}

#' Write out a completed (or partially-completed) golden-score run
#'
#' Writes `results.csv`, `metrics_permutations.csv`, `metrics_observed.csv`,
#' the faceted diagnostic histograms (see [plot_network_metric_hist_facet()]),
#' and `log_time.txt`, then deletes the per-condition temp files used for
#' resumability.
#'
#' @param results Results data frame accumulated so far.
#' @param all_metrics_df Long-format permutation-value rows accumulated so far.
#' @param all_obs_metrics_df Long-format observed-value rows accumulated so far.
#' @param temp_files Character vector of temp file paths to delete.
#' @param respath Results folder.
#' @param loop_start `Sys.time()` value from when the run started, for timing.
#' @param logs Character vector of log lines to append summary info to.
#'
#' @return A list: `results`, `logs`.
#' @export
finalize_golden_score_results <- function(results, all_metrics_df, all_obs_metrics_df, temp_files, respath, loop_start, logs) {
  readr::write_csv(results, paste0(respath, "/results.csv"))
  readr::write_csv(all_metrics_df, file.path(respath, "metrics_permutations.csv"))
  readr::write_csv(all_obs_metrics_df, file.path(respath, "metrics_observed.csv"))
  plot_network_metric_hist_facet(all_metrics_df, all_obs_metrics_df, respath)

  elapsed <- Sys.time() - loop_start
  print(paste0("Run time: ", elapsed))
  writeLines(as.character(elapsed), file.path(respath, "log_time.txt"))

  logs <- c(logs, paste0("Run time: ", elapsed), utils::capture.output(print(results)))

  for (f in temp_files) {
    if (file.exists(f)) file.remove(f)
  }

  list(results = results, logs = logs)
}

#' Reconstruct partial network-significance results from per-condition temp files
#'
#' Recovers `score_sig_network` values from `temp_network_<condition>_<perc_cutoff>.txt`
#' files left behind by an interrupted or crashed run (these are deleted by
#' [finalize_golden_score_results()] once a run completes normally).
#'
#' @param respath Results folder to look in.
#' @param conditions Character vector of condition/cell labels to check for.
#' @param perc_cutoffs Numeric vector of percentile cutoffs to check for.
#'
#' @return A data frame, one row per `condition` x `perc_cutoff` combination,
#'   with `score_sig_network` filled in where a temp file was found (`NA`
#'   otherwise).
#' @export
reconstruct_golden_score_from_temp <- function(respath, conditions, perc_cutoffs) {
  df <- expand.grid(condition = conditions, perc_cutoff = perc_cutoffs, stringsAsFactors = FALSE)
  df$score_sig_network <- NA_real_

  for (i in seq_len(nrow(df))) {
    f <- file.path(respath, paste0("temp_network_", df$condition[i], perc_suffix(df$perc_cutoff[i]), ".txt"))
    if (file.exists(f)) {
      val <- suppressWarnings(as.numeric(readLines(f, n = 1, warn = FALSE)))
      if (!is.na(val)) df$score_sig_network[i] <- val
    }
  }
  df
}

#' Find condition/perc_cutoff combinations with no temp file (never completed)
#'
#' @inheritParams reconstruct_golden_score_from_temp
#'
#' @return A data frame subset of `condition` x `perc_cutoff` combinations
#'   whose `temp_network_*.txt` file does not exist.
#' @export
find_missing_golden_score_combinations <- function(respath, conditions, perc_cutoffs) {
  df <- expand.grid(condition = conditions, perc_cutoff = perc_cutoffs, stringsAsFactors = FALSE)
  missing <- !file.exists(file.path(respath, paste0("temp_network_", df$condition, perc_suffix(df$perc_cutoff), ".txt")))
  df[missing, , drop = FALSE]
}
