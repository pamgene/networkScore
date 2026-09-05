test_that("clean_uka_to_kinograte_kinase selects and renames columns", {
  uka <- data.frame(
    `x.Sample` = "cond_A",
    `x.Kinase Name` = "KIN1",
    `x.Median Kinase Statistic` = 1.5,
    `x.Mean Specificity Score` = 2.0,
    check.names = FALSE
  )
  res <- clean_uka_to_kinograte_kinase(uka, cs = FALSE)
  expect_equal(colnames(res), c("Sample", "uniprotname", "LogFC", "fscore"))
})

test_that("clean_uka_to_kinograte_full flips sign when control is on the left of the contrast", {
  uka <- data.frame(
    `x.contrast` = c("Control vs Treated", "Treated2 vs Control"),
    `x.Kinase Name` = c("K1", "K2"),
    `x.Median Kinase Statistic` = c(1.0, 1.0),
    `x.Mean Specificity Score` = c(2.0, 2.0),
    check.names = FALSE
  )
  res <- clean_uka_to_kinograte_full(uka, spec_cutoff = 0, control = "Control", cs = FALSE)

  expect_equal(res$LogFC[res$cell_line == "Treated"], -1.0) # control on left -> sign flipped
  expect_equal(res$LogFC[res$cell_line == "Treated2"], 1.0) # control on right -> unchanged
})

test_that("clean_sens_to_kinograte computes fold-change vs control and can restrict to best_drug_per_target", {
  sens <- data.frame(
    CELL_LINE_NAME = c("Control", "Control", "A", "A"),
    TARGET_1 = c("T1", "T2", "T1", "T2"),
    LN_IC50 = c(1.0, 2.0, 3.0, 1.5),
    DRUG_ID = c("D1", "D2", "D1", "D2")
  )
  res <- clean_sens_to_kinograte(sens, control = "Control")
  expect_equal(res$LogFC[res$uniprotname == "T1"], 2.0) # 3.0 - 1.0
  expect_equal(res$LogFC[res$uniprotname == "T2"], -0.5) # 1.5 - 2.0
  expect_false("Control" %in% res$cell_line)

  res_restricted <- clean_sens_to_kinograte(sens, control = "Control", best_drug_per_target = data.frame(DRUG_ID = "D1"))
  expect_equal(unique(res_restricted$uniprotname), "T1")
})
