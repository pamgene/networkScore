#' Clean a raw UKA table for the paired (golden) analysis, keyed by cell line/contrast
#'
#' Handles both "X vs control" and "control vs X" contrast naming, flipping
#' the statistic's sign when control is on the left. Unlike
#' `networkGen::clean_uka_to_kinograte1()` (which this is a sibling of), the
#' output is keyed `cell_line` rather than `Sgroup_contrast`, matching what
#' [make_golden_score_full()] expects.
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param spec_cutoff Minimum specificity-score to keep a row.
#' @param control Name of the control condition as it appears in `contrast`.
#' @param cs `TRUE`/`FALSE` to force per-comparison (csUKA) vs. mean/median
#'   columns; `NULL` (default) auto-detects via [networkGen::detect_csuka()].
#'
#' @return Data frame with columns `cell_line`, `comparison` (the full,
#'   untruncated "X vs control"/"control vs X" label -- `cell_line` itself
#'   has `control` stripped out, since it's used as a join key against
#'   sensitivity data's own cell-line identifier, not shown to the user),
#'   `uniprotname`, `LogFC`, `fscore`.
#' @export
clean_uka_to_kinograte_full <- function(uka, spec_cutoff, control, cs = NULL) {
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

#' Clean a raw UKA table for the kinase-only analysis, keyed by condition
#'
#' The kinase-only counterpart of `networkGen::clean_uka_to_kinograte()`:
#' same reduction of a raw Tercen UKA export to the shape the grid builder
#' consumes, differing only in which raw column identifies each condition
#' (`condition_col`).
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param cs `TRUE`/`FALSE` to force per-comparison (csUKA) vs. mean/median
#'   columns; `NULL` (default) auto-detects via [networkGen::detect_csuka()].
#' @param condition_col Name of the raw column identifying each condition/
#'   comparison, carried through unchanged to the output. This is the same
#'   role `Sgroup_contrast` plays in `networkGen::clean_uka_to_kinograte()`;
#'   different Tercen exports name it differently (`"Sgroup_contrast"`,
#'   `"Sample"`, ...). Default `"Sgroup_contrast"`. Must match the
#'   `condition_col` passed to `networkGen::build_network_grid()`.
#'
#' @return Data frame with columns `<condition_col>`, `uniprotname`,
#'   `LogFC`, `fscore`.
#' @export
clean_uka_to_kinograte_kinase <- function(uka, cs = NULL, condition_col = "Sgroup_contrast") {
  if (is.null(cs)) cs <- networkGen::detect_csuka(uka)
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  stat_col <- if (cs) "Kinase Statistic" else "Median Kinase Statistic"

  uka %>%
    networkGen::clean_tercen_columns() %>%
    dplyr::select(dplyr::all_of(condition_col), "Kinase Name", dplyr::all_of(stat_col), dplyr::all_of(finalscore_col)) %>%
    dplyr::rename(
      "uniprotname" = "Kinase Name", "LogFC" = dplyr::all_of(stat_col),
      "fscore" = dplyr::all_of(finalscore_col)
    ) %>%
    dplyr::distinct()
}

#' Clean sensitivity data for the paired (golden) analysis
#'
#' Sibling of `networkGen::clean_sens_to_kinograte()`; this variant renames
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
clean_sens_to_kinograte <- function(sens, control, zscore = FALSE, del_cell = NULL, best_drug_per_target = NULL) {
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

