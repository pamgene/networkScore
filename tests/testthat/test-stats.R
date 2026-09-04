fake_result <- function(edges_df, missing_nodes = NULL) {
  g <- igraph::graph_from_data_frame(edges_df, directed = FALSE)
  structure(list(network = g, missing_nodes = missing_nodes), class = "networkGen_result")
}

test_that("compute_network_stats returns the four expected fields, no modularity/assortativity/core", {
  result <- fake_result(data.frame(from = c("A", "B"), to = c("B", "C")))
  stats <- compute_network_stats(result)

  expect_named(stats, c("rel_med_path", "rel_med_path_inv", "density", "clustering"))
  expect_true(is.numeric(stats$rel_med_path))
  expect_true(is.numeric(stats$density) && stats$density >= 0 && stats$density <= 1)
})

test_that("rel_med_path_inv is the inverse-transformed version -- tighter networks score higher", {
  tight <- fake_result(data.frame(from = c("A", "A", "A"), to = c("B", "C", "D")))
  loose <- fake_result(data.frame(from = c("A", "B", "C"), to = c("B", "C", "D"))) # a path, more spread out

  tight_stats <- compute_network_stats(tight)
  loose_stats <- compute_network_stats(loose)

  expect_true(tight_stats$rel_med_path_inv >= loose_stats$rel_med_path_inv)
})

test_that("compute_network_stats folds missing_nodes into n_nodes and pads distances at the observed max", {
  # A-B: distances(mode="all") as a vector is c(0,1,1,0), median 0.5, n_nodes 2 -> rel_med_path 0.25
  result_no_missing <- fake_result(data.frame(from = "A", to = "B"))
  s1 <- compute_network_stats(result_no_missing)
  expect_equal(s1$rel_med_path, 0.25)

  # + 2 missing nodes: max_dist (1) appended twice -> c(0,1,1,0,1,1), median 1, n_nodes 4 -> 0.25
  result_with_missing <- fake_result(data.frame(from = "A", to = "B"), missing_nodes = data.frame(name = c("X", "Y")))
  s2 <- compute_network_stats(result_with_missing)
  expect_equal(s2$rel_med_path, 0.25)

  # + 5 missing nodes: median shifts further toward max_dist while n_nodes grows faster -> net decrease
  result_many_missing <- fake_result(data.frame(from = "A", to = "B"), missing_nodes = data.frame(name = paste0("X", 1:5)))
  s3 <- compute_network_stats(result_many_missing)
  expect_true(s3$rel_med_path < s1$rel_med_path)
})
