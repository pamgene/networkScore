#' @importFrom rlang .data
#' @importFrom dplyr %>%
NULL

# Bare column-name symbols passed positionally into networkGen's tidy-eval
# helpers (percentile_score_fast()/percentile_score_noabs()'s `symbol`/`metric`
# args) -- not real global variables, but R CMD check's static analysis
# can't tell the difference.
utils::globalVariables(c("uniprotname", "LogFC"))
