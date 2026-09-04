#' Compute network-topology statistics for a scored network
#'
#' The statistics used to compare an observed network against its permuted
#' null distribution. This is the reconciled replacement for
#' `Network_generation`'s `compute_network_metrics_kinase()` (the only one of
#' the two original statistics functions that was ever actually implemented
#' -- the paired-data variant was commented out and never finished) --
#' ported unchanged in substance, used for both the paired and kinase-only
#' scoring paths (both build a `networkGen_result`, and this function only
#' looks at its `network`/`missing_nodes`, so one implementation serves
#' both).
#'
#' `modularity`/`assortativity` (present as placeholder fields elsewhere in
#' the pre-migration code, never actually computed) and the `_core`/`_core_inv`
#' fields the old permutation code read but nothing computed are deliberately
#' not reproduced here -- see `docs/adr/` at the `DevOpti` root.
#'
#' @param network_result A `"networkGen_result"` object (from
#'   `networkGen::generate_paired_network()` / `generate_kinase_network()`).
#' @param a Smoothing constant in the inverse shortest-path transform
#'   (`1 / (a + distance)`). Default 1.
#'
#' @return A list: `rel_med_path` (median shortest-path length across all
#'   node pairs, relative to node count -- lower means a tighter network),
#'   `rel_med_path_inv` (same, inverse-transformed so higher means tighter),
#'   `density`, `clustering` (mean local transitivity).
#' @export
compute_network_stats <- function(network_result, a = 1) {
  g <- network_result$network
  missing_nodes <- network_result$missing_nodes

  shortest_path_stat <- function(inverse = FALSE) {
    all_dist <- as.vector(igraph::distances(g, mode = "all"))
    all_dist[is.infinite(all_dist)] <- NA

    if (!is.null(missing_nodes) && nrow(missing_nodes) > 0) {
      n_nodes <- igraph::vcount(g) + nrow(missing_nodes)
      max_dist <- max(all_dist, na.rm = TRUE)
      all_dist <- c(all_dist, rep(max_dist, nrow(missing_nodes)))
    } else {
      n_nodes <- igraph::vcount(g)
    }

    if (inverse) {
      all_dist <- 1 / (a + all_dist)
    }

    stats::median(all_dist, na.rm = TRUE) / n_nodes
  }

  list(
    rel_med_path = shortest_path_stat(inverse = FALSE),
    rel_med_path_inv = shortest_path_stat(inverse = TRUE),
    density = igraph::edge_density(g, loops = FALSE),
    clustering = mean(igraph::transitivity(g, type = "local"), na.rm = TRUE)
  )
}
