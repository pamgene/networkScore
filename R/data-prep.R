#' Reshape a raw UKA table for the paired (golden) analysis, keyed by cell line
#'
#' The paired counterpart of `networkGen::prep_uka()`. Given a fixed
#' `control`, it keeps only the comparisons that involve that control
#' (either "X vs control" or "control vs X"), flips the kinase statistic's
#' sign when the control is on the left so `LogFC` always reads "treatment
#' relative to control", and derives `cell_line` -- the non-control side of
#' each comparison -- as the join key against the sensitivity data. Any
#' comparison not involving `control` is dropped.
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param spec_cutoff Minimum specificity-score to keep a row.
#' @param control Name of the control condition as it appears in `contrast`
#'   (e.g. `"DMSO"`). Set this to whichever side of every comparison is the
#'   reference; the other side becomes `cell_line`.
#' @param cs `TRUE`/`FALSE` to force per-comparison (csUKA) vs. mean/median
#'   columns; `NULL` (default) auto-detects via [networkGen::detect_csuka()].
#'
#' @return Data frame with columns `cell_line`, `comparison` (the full,
#'   untruncated "X vs control"/"control vs X" label -- `cell_line` itself
#'   has `control` stripped out, since it's used as a join key against
#'   sensitivity data's own cell-line identifier, not shown to the user),
#'   `uniprotname`, `LogFC`, `fscore`.
#' @export
prep_uka_paired <- function(uka, spec_cutoff, control, cs = NULL) {
  if (is.null(cs)) cs <- networkGen::detect_csuka(uka)
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  stat_col <- if (cs) "Kinase Statistic" else "Median Kinase Statistic"

  uka %>%
    networkGen::clean_tercen_columns() %>%
    dplyr::filter(.data[[finalscore_col]] > spec_cutoff) %>%
    dplyr::filter(
      stringr::str_detect(.data$contrast, paste0(" vs ", control)) |
        stringr::str_detect(.data$contrast, paste0(control, " vs "))
    ) %>%
    dplyr::mutate(contrast = gsub(.data$contrast, pattern = "-", replacement = "")) %>%
    dplyr::mutate(
      control_left = stringr::str_detect(.data$contrast, paste0("^", control, " vs ")),
      comparison = .data$contrast,
      cell_line = ifelse(.data$control_left,
        stringr::str_replace(.data$contrast, paste0(control, " vs "), ""),
        stringr::str_replace(.data$contrast, paste0(" vs ", control), "")
      ),
      MKS = ifelse(.data$control_left, -.data[[stat_col]], .data[[stat_col]])
    ) %>%
    dplyr::select("cell_line", "comparison", "Kinase Name", "MKS", dplyr::all_of(finalscore_col)) %>%
    dplyr::rename("uniprotname" = "Kinase Name", "LogFC" = "MKS", "fscore" = dplyr::all_of(finalscore_col)) %>%
    dplyr::mutate(cell_line = gsub(.data$cell_line, pattern = "-", replacement = "")) %>%
    dplyr::distinct()
}

#' Reshape sensitivity data for the paired (golden) analysis
#'
#' Sibling of `networkGen::prep_sens()`; this variant renames
#' `CELL_LINE_NAME` (rather than expecting a pre-named `cell_line` column)
#' and can additionally restrict to a preferred drug per target.
#'
#' @param sens Raw sensitivity data frame with columns `CELL_LINE_NAME`,
#'   `TARGET_1`, `LN_IC50` (and `DRUG_ID` if `best_drug_per_target` is used).
#' @param control If given, `LogFC` is fold-change vs. this `cell_line`. If
#'   `NULL` and `zscore = FALSE`, `sens` is returned as-is. If `NULL` and
#'   `zscore = TRUE`, a per-target z-score across cell lines is computed.
#' @param zscore See `control`.
#' @param del_cell Optional character vector of `cell_line` values to drop.
#' @param best_drug_per_target Optional data frame with a `DRUG_ID` column;
#'   if given, `sens` is restricted to rows whose `DRUG_ID` is in it.
#'
#' @return Data frame with columns `cell_line`, `uniprotname`, `LogFC`.
#' @export
prep_sens <- function(sens, control, zscore = FALSE, del_cell = NULL, best_drug_per_target = NULL) {
  sens_filt <- sens %>% dplyr::rename(cell_line = "CELL_LINE_NAME", uniprotname = "TARGET_1")

  if (!is.null(best_drug_per_target)) {
    sens_filt <- sens_filt %>% dplyr::filter(.data$DRUG_ID %in% best_drug_per_target$DRUG_ID)
  }

  if (!is.null(control) && !zscore) {
    control_df <- sens_filt %>% dplyr::filter(.data$cell_line == control)
    sens_clean <- sens_filt %>%
      dplyr::left_join(control_df, by = c("uniprotname"), suffix = c("", ".control")) %>%
      dplyr::mutate(LogFC = .data$LN_IC50 - .data$LN_IC50.control) %>%
      dplyr::filter(.data$cell_line != control, !is.na(.data$LogFC)) %>%
      dplyr::select("cell_line", "uniprotname", "LogFC")
  } else if (is.null(control) && !zscore) {
    sens_clean <- sens_filt
  } else if (is.null(control) && zscore) {
    sens_clean <- sens_filt %>%
      tidyr::pivot_wider(id_cols = "cell_line", names_from = "uniprotname", values_from = "LN_IC50") %>%
      tibble::column_to_rownames("cell_line") %>%
      as.matrix() %>%
      scale() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell_line") %>%
      tidyr::pivot_longer(cols = -"cell_line", names_to = "uniprotname", values_to = "LogFC")
  }

  if (!is.null(del_cell)) {
    sens_clean <- sens_clean %>% dplyr::filter(!.data$cell_line %in% del_cell)
  }

  sens_clean
}

