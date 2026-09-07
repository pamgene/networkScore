#' @importFrom rlang .data
#' @importFrom dplyr %>%
NULL

# Bare column-name symbols passed positionally into networkGen's tidy-eval
# helpers (percentile_score_fast()/percentile_score_noabs()'s `symbol`/`metric`
# args) -- not real global variables, but R CMD check's static analysis
# can't tell the difference.
utils::globalVariables(c("uniprotname", "LogFC"))

#' Build the `_<perc_cutoff>` temp-file suffix, omitted when there's no filter
#'
#' `perc_cutoff = 0` means "no percentile filtering" -- encoding it in every
#' temp-file name would just be noise, so it's dropped entirely in that case.
#' Used both when writing temp files and when looking them back up (see
#' [reconstruct_golden_score_from_temp()]/[find_missing_golden_score_combinations()])
#' -- both sides must apply the same rule or they'll disagree on the filename.
#'
#' @param perc_cutoff The percentile cutoff used for this run.
#'
#' Vectorized: safe to use on a single value or a vector of them (see
#' [find_missing_golden_score_combinations()]).
#'
#' @return `""` where `perc_cutoff` is `0`; `"_<perc_cutoff>"` elsewhere.
#' @keywords internal
perc_suffix <- function(perc_cutoff) {
  ifelse(perc_cutoff == 0, "", paste0("_", perc_cutoff))
}
