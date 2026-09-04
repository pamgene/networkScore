test_that("pOverlapScore computes a floored permutation significance", {
  set.seed(1)
  uka_cell_all <- data.frame(uniprotname = paste0("K", 1:20), LogFC = c(rep(3, 5), rnorm(15)), fscore = 2)
  sens_filt <- data.frame(name = paste0("K", 1:5)) # perfectly overlaps the top 5 kinases by |LogFC|

  score_obs <- networkGen::overlap_uka_sens(paste0("K", 1:5), sens_filt$name)
  sig <- pOverlapScore(
    score_obs = score_obs, uka_cell_all = uka_cell_all, sens_filt = sens_filt, cell = "test",
    nsample = 5, nPerms = 20, spec_cutoff = 0, rank_uka_abs = TRUE, perc_cutoff = 0.7
  )

  expect_true(sig >= 1 / 20 && sig <= 1)
})

test_that("compute_golden_overlap_scores writes temp files and returns both individual and family-wise scores", {
  set.seed(1)
  respath <- file.path(tempdir(), "overlap_test")
  dir.create(respath, showWarnings = FALSE)

  uka_filt <- data.frame(name = c("K1", "K2"))
  sens_filt <- data.frame(name = c("K1", "K3"))
  uka_cell_all <- data.frame(uniprotname = c("K1", "K2", "K3", "K4"), LogFC = c(3, 2, 1, 0.1), fscore = 2)
  uka_fam <- data.frame(Kinase_Name = c("K1", "K2", "K3", "K4"), Kinase_family = c("FamA", "FamA", "FamB", "FamB"))

  res <- compute_golden_overlap_scores(
    uka_filt = uka_filt, sens_filt = sens_filt, uka_cell_all = uka_cell_all, uka_fam = uka_fam,
    cell = "test_cell", perc_cutoff = 0.5, nperms_overlap = 10, respath = respath,
    spec_cutoff = 0, rank_uka_abs = TRUE
  )

  expect_named(res$vals, c("obs_overlap", "score_sig_overlap", "obs_overlap_fam", "score_sig_overlap_fam"))
  expect_equal(res$vals$obs_overlap, 0.5) # 1 of 2 sens names overlap
  expect_length(res$temp_files, 2)
  expect_true(all(file.exists(res$temp_files)))
})
