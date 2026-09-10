#' Extend networkGen's folder-naming rule with networkScore-only parameters
#'
#' Composes [networkGen::build_param_folder()] (the `networkGen`-native
#' components: `uka`/`sens`/`ppi_network`/`spec_cutoff`/`perc_cutoff`/
#' `rank_uka_abs`/`b`/`w`/`cs`/`art_nodes`) with the scoring-only components
#' this package's own naming needs -- `nperms_network`, `relative_to`.
#' `networkGen` deliberately has no knowledge of these (it has no scoring
#' concept at all -- see the suite's `CONTEXT-MAP.md`), so this composes on
#' top of its naming rule rather than duplicating it.
#'
#' @param params A parameters list from [networkGen::capture_params()].
#'
#' @return The folder name, as a string.
#' @export
build_score_param_folder <- function(params) {
  base <- networkGen::build_param_folder(params)
  extra <- c(
    if ("nperms_network" %in% names(params)) paste0("nperms", params$nperms_network),
    if ("relative_to" %in% names(params)) paste0("rel2", params$relative_to)
  )
  if (length(extra) == 0) base else paste(c(base, extra), collapse = "_")
}

#' Capture a scoring run's parameters, build its output folder, and write `params.csv`
#'
#' Composes [networkGen::capture_params()] and [build_score_param_folder()],
#' then the filesystem side effects -- mirrors
#' [networkGen::prepare_run_params()] but uses this package's extended
#' naming rule instead of `networkGen`'s alone. Used by
#' [make_golden_score_kinase()] / [make_golden_score_full()].
#'
#' @param ... Passed to [networkGen::capture_params()]. Must include
#'   `respath` (the base output directory this run's parameter-encoded
#'   folder is created under) and whatever [build_score_param_folder()]
#'   needs.
#'
#' @return The evaluated parameters as a named list, with `respath` replaced
#'   by the created parameter-encoded folder path.
#' @export
prepare_score_run_params <- function(...) {
  params <- networkGen::capture_params(...)
  folder_name <- build_score_param_folder(params)
  param_folder <- file.path(params$respath, folder_name)
  print(paste0("param folder: ", param_folder))
  if (!dir.exists(param_folder)) dir.create(param_folder, recursive = TRUE)
  params$respath <- param_folder
  networkGen::save_params(params, respath = param_folder)
  params
}
