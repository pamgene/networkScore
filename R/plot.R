#' Faceted permutation-null-distribution diagnostic plot
#'
#' One PNG per `metric` x `perc_cutoff` combination, faceted by condition:
#' each subplot is a histogram of that condition's permutation values for
#' that metric, with the observed value marked by a vertical line. Lets you
#' scan many conditions at once for which ones have their observed value out
#' in the tail (significant) vs. buried in the null distribution.
#'
#' This is a scoring diagnostic, not network visualization -- it belongs
#' here in `networkScore`, not in `networkPlot`. Migrated from
#' `Network_generation`'s `R/golden_score_helper.R`, which already
#' generalized the original per-project (kinase-only vs. paired) split by
#' detecting whether the input has a `cell` or `condition` grouping column;
#' kept as-is here.
#'
#' @param all_metrics_df Long-format permutation-value rows (`metric`,
#'   `perc_cutoff`, `value`, and a `cell` or `condition` grouping column).
#' @param all_obs_df Long-format observed-value rows, same shape.
#' @param res.path Folder to write `Network_metrics_plots/<metric>_cutoff<cutoff>.png` into.
#'
#' @return Invisibly, the list of `ggplot` objects produced (also saved to disk).
#' @export
plot_network_metric_hist_facet <- function(all_metrics_df, all_obs_df, res.path) {
  dir.create(file.path(res.path, "Network_metrics_plots"), showWarnings = FALSE, recursive = TRUE)

  if (nrow(all_metrics_df) == 0) {
    return(invisible(list()))
  }

  make_plot <- function(metric_i, cutoff_i) {
    df1_sub <- all_metrics_df %>%
      dplyr::filter(.data$metric == metric_i, .data$perc_cutoff == cutoff_i, is.finite(.data$value), !is.na(.data$value))
    df2_sub <- all_obs_df %>%
      dplyr::filter(.data$metric == metric_i, .data$perc_cutoff == cutoff_i, is.finite(.data$value), !is.na(.data$value))

    fname <- paste0(res.path, "/Network_metrics_plots/", metric_i, "_cutoff", cutoff_i, ".png")

    if (nrow(df1_sub) == 0 || nrow(df2_sub) == 0) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data available", size = 6) +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0("Metric: ", metric_i, " | Cutoff: ", cutoff_i))
      ggplot2::ggsave(fname, p, width = 7, height = 5, dpi = 300)
      return(p)
    }

    grouping_col <- if ("cell" %in% colnames(df1_sub)) "cell" else "condition"

    p <- ggplot2::ggplot(df1_sub, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "black") +
      ggplot2::geom_vline(data = df2_sub, ggplot2::aes(xintercept = .data$value), color = "red", linetype = "dashed") +
      ggplot2::facet_wrap(stats::as.formula(paste("~", grouping_col)), scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "top",
        strip.background = ggplot2::element_rect(fill = "grey90"),
        strip.text = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(title = paste0("Metric: ", metric_i, " | Cutoff: ", cutoff_i), x = "Permutation values", y = "Count")

    ggplot2::ggsave(fname, p, width = 7, height = 5, dpi = 300)
    p
  }

  combos <- all_metrics_df %>% dplyr::distinct(.data$metric, .data$perc_cutoff)
  invisible(purrr::pmap(list(combos$metric, combos$perc_cutoff), make_plot))
}
