#' Permutation significance for the (family-wise) kinase/sensitivity overlap
#'
#' Repeatedly shuffles kinase-activity values across kinases, re-selects top
#' hits, and recomputes the overlap with the sensitivity hits, to build a
#' null distribution for `score_obs`. Cheap relative to network-based
#' scoring (no PCSF calls), so this stays sequential rather than going
#' through `networkGen::generate_networks_batch()`.
#'
#' @param score_obs The observed overlap score (from `networkGen::overlap_uka_sens()`).
#' @param uka_cell_all All kinase-activity rows for this condition (unfiltered).
#' @param sens_filt Top-hit sensitivity data frame for this condition.
#' @param cell Condition/cell-line label, for logging.
#' @param nsample Unused, kept for call-site compatibility.
#' @param nPerms Number of permutations. Default 1000.
#' @param spec_cutoff,rank_uka_abs,perc_cutoff Passed to `networkGen::uka_top()`.
#' @param uka_fam Kinase-family lookup table (`Kinase_Name`, `Kinase_family`
#'   columns), required when `family_mode = TRUE`.
#' @param family_mode If `TRUE`, compute the overlap at kinase-family level
#'   instead of individual-kinase level.
#'
#' @return The significance score: the fraction of permutations whose
#'   overlap was at least as large as `score_obs` (floored at `1 / nPerms`).
#' @export
pOverlapScore <- function(score_obs, uka_cell_all, sens_filt, cell, nsample, nPerms = 1000,
                           spec_cutoff, rank_uka_abs, perc_cutoff, uka_fam = NULL, family_mode = FALSE) {
  score_perm <- rep(NA, nPerms)

  for (i in seq_len(nPerms)) {
    uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
    uka_filt_random <- uka_cell_all %>% networkGen::uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)

    if (!family_mode) {
      score_perm[i] <- networkGen::overlap_uka_sens(uka_filt_random$name, sens_filt$name)
    } else {
      uka_fam_rand <- uka_fam %>% dplyr::filter(.data$Kinase_Name %in% uka_filt_random$name) %>% dplyr::pull(.data$Kinase_family) %>% unique()
      sens_fam <- uka_fam %>% dplyr::filter(.data$Kinase_Name %in% sens_filt$name) %>% dplyr::pull(.data$Kinase_family) %>% unique()
      score_perm[i] <- networkGen::overlap_uka_sens(uka_fam_rand, sens_fam)
    }
  }

  overlap_sig <- max(sum(score_perm >= score_obs) / nPerms, 1 / nPerms)
  print(paste0("Cell: ", cell, ", overlap score", ifelse(family_mode, " family: ", ": "), overlap_sig))
  overlap_sig
}

#' Observed overlap score and its permutation significance, individual + family-wise
#'
#' @param uka_filt Top-hit kinase data frame for this condition.
#' @param sens_filt Top-hit sensitivity data frame for this condition.
#' @param uka_cell_all All kinase-activity rows for this condition (unfiltered).
#' @param uka_fam Kinase-family lookup table (`Kinase_Name`, `Kinase_family`).
#' @param cell Condition/cell-line label.
#' @param perc_cutoff,spec_cutoff,rank_uka_abs Passed through to [pOverlapScore()].
#' @param nperms_overlap Number of permutations for both overlap scores.
#' @param respath Folder to write temp progress files into (deleted by
#'   [finalize_golden_score_results()] once the whole run completes; lets a
#'   crashed/interrupted run be resumed).
#'
#' @return A list: `vals` (named list with `obs_overlap`, `score_sig_overlap`,
#'   `obs_overlap_fam`, `score_sig_overlap_fam`) and `temp_files` (paths written).
#' @export
compute_golden_overlap_scores <- function(uka_filt, sens_filt, uka_cell_all, uka_fam, cell,
                                           perc_cutoff, nperms_overlap, respath, spec_cutoff, rank_uka_abs) {
  vals <- list()

  score_obs <- networkGen::overlap_uka_sens(uka_filt$name, sens_filt$name)
  vals$obs_overlap <- round(score_obs, 2)

  overlap_sig <- pOverlapScore(
    score_obs = score_obs, uka_cell_all = uka_cell_all, sens_filt = sens_filt, cell = cell,
    nsample = nrow(uka_filt), nPerms = nperms_overlap, spec_cutoff = spec_cutoff,
    rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff
  )
  vals$score_sig_overlap <- overlap_sig

  temp_overlap_file <- file.path(respath, paste0("temp_overlap_", cell, "_", perc_cutoff, ".txt"))
  writeLines(as.character(overlap_sig), temp_overlap_file)
  temp_files <- c(temp_overlap_file)

  uka_filt_fam <- uka_fam %>% dplyr::filter(.data$Kinase_Name %in% uka_filt$name) %>% dplyr::pull(.data$Kinase_family) %>% unique()
  sens_filt_fam <- uka_fam %>% dplyr::filter(.data$Kinase_Name %in% sens_filt$name) %>% dplyr::pull(.data$Kinase_family) %>% unique()
  score_obs_fam <- networkGen::overlap_uka_sens(uka_filt_fam, sens_filt_fam)
  vals$obs_overlap_fam <- round(score_obs_fam, 2)

  overlap_sig_fam <- pOverlapScore(
    score_obs = score_obs_fam, uka_cell_all = uka_cell_all, sens_filt = sens_filt, cell = cell,
    nsample = nrow(uka_filt), nPerms = nperms_overlap, spec_cutoff = spec_cutoff,
    rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff, uka_fam = uka_fam, family_mode = TRUE
  )
  vals$score_sig_overlap_fam <- overlap_sig_fam

  temp_overlap_fam_file <- file.path(respath, paste0("temp_overlap_fam_", cell, "_", perc_cutoff, ".txt"))
  writeLines(as.character(overlap_sig_fam), temp_overlap_fam_file)
  temp_files <- c(temp_files, temp_overlap_fam_file)

  list(vals = vals, temp_files = temp_files)
}
