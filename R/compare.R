#' Compare network-topology significance across several result folders
#'
#' Reads every `results.csv` found (recursively -- so a folder can be a
#' single grid combination, or a whole multi-combination run) under each
#' given result folder, and plots `score_sig_network`/`score_sig_network_inv`
#' (by default) as boxplots per `Comparison`, one facet per folder, so the
#' *same* comparisons across *different* grid combinations -- different runs,
#' different `respath`s -- can be compared side by side. Within one facet,
#' each box's spread comes from however many grid combinations exist inside
#' that folder (e.g. several `b`/`spec_cutoff` subfolders under one base
#' `respath`) -- a single-combination folder gives a one-point "box".
#'
#' @param respaths A named character vector (or named list) of result
#'   folder paths -- names are the aliases shown as facet titles, values are
#'   the folder paths to search (each searched recursively for
#'   `results.csv` files, so either a single grid-combination folder or a
#'   base `respath` covering many combinations works). Facets appear in the
#'   order given here, not alphabetically.
#' @param metrics Character vector of score columns to compare, shown as
#'   different colors within each box group. Default `score_sig_network`/
#'   `score_sig_network_inv`.
#' @param sig_threshold Numeric, or `NULL` to omit. Draws a dashed
#'   horizontal reference line at this value (default `0.05`, the
#'   conventional significance cutoff for these permutation-test p-values).
#'
#' @return A `ggplot` object (not saved to disk -- call [ggplot2::ggsave()]
#'   on the result if you want a file).
#'
#' @details
#' Only `Comparison` values present in *every* listed folder are plotted --
#' the point of this plot is comparing the same comparisons across grid
#' combinations, so a comparison missing from even one folder would leave a
#' misleading gap in that folder's facet. If any folder has comparisons the
#' others don't, those are silently dropped and a message reports what was
#' excluded and from where. If no `Comparison` is common to every folder at
#' all, this is refused with `stop()` rather than silently producing an
#' empty plot -- almost always an accidental folder selection (e.g. two
#' unrelated runs), not an intended empty result.
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
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::facet_wrap(~alias) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(x = "Comparison", y = "score (p-value)", fill = "Metric")

  if (!is.null(sig_threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = sig_threshold, linetype = "dashed")
  }

  p
}
