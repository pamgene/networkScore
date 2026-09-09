#' Compare network-topology significance across several result folders
#'
#' Reads the `results.csv` under each given result folder and plots
#' `score_sig_network`/`score_sig_network_inv` (by default) as bars per
#' `Comparison`, one facet per folder -- for comparing the *same*
#' comparisons across *different* parameter combinations or entirely
#' separate runs, one folder per combination (matching how
#' [networkGen::run_network_grid()]/[make_golden_score_kinase()]/
#' [make_golden_score_full()] lay out output: one folder per distinct
#' combination).
#'
#' @param respaths A named character vector (or named list) of result
#'   folder paths -- names are the aliases shown as facet titles, values are
#'   the folder paths to read (each searched recursively for `results.csv`
#'   files, so pointing at a base `respath` covering several combinations
#'   also works, *provided* it doesn't produce more than one row for the
#'   same `Comparison` -- see Details). Facets appear in the order given
#'   here, not alphabetically.
#' @param metrics Character vector of score columns to compare, shown as
#'   different colors within each bar group. Default `score_sig_network`/
#'   `score_sig_network_inv`.
#' @param sig_threshold Numeric, or `NULL` to omit. Draws a dashed
#'   horizontal reference line at this value (default `0.05`, the
#'   conventional significance cutoff for these permutation-test p-values).
#'
#' @return A `ggplot` object (not saved to disk -- call [ggplot2::ggsave()]
#'   on the result if you want a file). The object has no size of its own --
#'   with several folders and/or many comparisons, the default plot device
#'   size (e.g. a knitr chunk's default `fig.width`/`fig.height`) is usually
#'   too small and crowds the facets/labels; size explicitly for the number
#'   of folders and comparisons involved, e.g. `fig.width = 4 * n_folders`
#'   in a chunk, or `ggplot2::ggsave(..., width = 4 * n_folders, height = 6)`.
#'
#' @details
#' Each folder is expected to represent exactly *one* parameter combination
#' -- one row per `Comparison` -- so each bar has one real value behind it,
#' not an aggregate. If a folder's `results.csv` (or several found under it)
#' has more than one row for the same `Comparison`, this is refused with a
#' clear `stop()` naming the folder and the repeated comparison, rather than
#' silently plotting a bar that doesn't represent a single combination.
#'
#' Only `Comparison` values present in *every* listed folder are plotted --
#' the point of this plot is comparing the same comparisons across
#' combinations, so a comparison missing from even one folder would leave a
#' misleading gap in that folder's facet. If any folder has comparisons the
#' others don't, those are dropped and a message reports what was excluded
#' and from where. If no `Comparison` is common to every folder at all, this
#' is refused with `stop()` rather than silently producing an empty plot --
#' almost always an accidental folder selection (e.g. two unrelated runs),
#' not an intended empty result.
#'
#' @export
plot_score_comparison <- function(respaths, metrics = c("score_sig_network", "score_sig_network_inv"), sig_threshold = 0.05) {
  if (is.null(names(respaths)) || any(names(respaths) == "")) {
    stop(
      "respaths must be a named vector/list -- each name is the alias shown ",
      "as that folder's facet title, e.g. c(baseline = \"path/to/run1\", ",
      "with_KL = \"path/to/run2\").",
      call. = FALSE
    )
  }

  per_alias <- purrr::imap(respaths, function(path, alias) {
    files <- list.files(path, pattern = "^results\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(files) == 0) {
      stop("No results.csv found under the '", alias, "' folder: ", path, call. = FALSE)
    }

    df <- dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))
    missing_metrics <- setdiff(metrics, colnames(df))
    if (length(missing_metrics) > 0) {
      stop(
        "The '", alias, "' folder's results.csv has no ",
        paste(missing_metrics, collapse = "/"), " column(s) -- found: ",
        paste(colnames(df), collapse = ", "), ".",
        call. = FALSE
      )
    }

    dupes <- unique(df$Comparison[duplicated(df$Comparison)])
    if (length(dupes) > 0) {
      stop(
        "The '", alias, "' folder has more than one row for comparison ",
        paste(sprintf('"%s"', dupes), collapse = ", "), " (", length(files), " results.csv file(s) found under ",
        path, "). plot_score_comparison() expects one folder per parameter ",
        "combination -- point to that specific combination's folder instead ",
        "of one covering more than one.",
        call. = FALSE
      )
    }

    df$alias <- alias
    df
  })

  comparisons_by_alias <- purrr::map(per_alias, ~ unique(.x$Comparison))
  common <- Reduce(intersect, comparisons_by_alias)

  if (length(common) == 0) {
    detail <- paste(
      sprintf("  %s: %s", names(comparisons_by_alias), purrr::map_chr(comparisons_by_alias, paste, collapse = ", ")),
      collapse = "\n"
    )
    stop(
      "No Comparison is common to every selected result folder -- nothing to compare. ",
      "Comparisons found per folder:\n", detail,
      call. = FALSE
    )
  }

  excluded <- purrr::imap(comparisons_by_alias, ~ setdiff(.x, common))
  excluded <- excluded[lengths(excluded) > 0]
  if (length(excluded) > 0) {
    message(
      "Dropping comparisons not present in every selected folder (kept only ",
      "those common to all): ",
      paste(sprintf("%s has extra %s", names(excluded), purrr::map_chr(excluded, paste, collapse = "/")), collapse = "; ")
    )
  }

  combined <- dplyr::bind_rows(per_alias) %>%
    dplyr::filter(.data$Comparison %in% common) %>%
    dplyr::mutate(alias = factor(.data$alias, levels = names(respaths))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(metrics), names_to = "metric", values_to = "value")

  p <- ggplot2::ggplot(combined, ggplot2::aes(x = .data$Comparison, y = .data$value, fill = .data$metric)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, color = "black") +
    ggplot2::facet_wrap(~alias, labeller = ggplot2::label_wrap_gen(width = 25)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    ggplot2::labs(x = "Comparison", y = "score (p-value)", fill = "Metric")

  if (!is.null(sig_threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = sig_threshold, linetype = "dashed")
  }

  p
}
