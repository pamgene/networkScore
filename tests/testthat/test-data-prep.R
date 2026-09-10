test_that("prep_uka_paired auto-detects csUKA columns when cs is NULL (default)", {
  cs_full <- data.frame(
    `x.contrast` = "Treated vs Control", `x.Kinase Name` = "K1",
    `x.Kinase Statistic` = 1.0, `x.Specificity Score` = 2.0,
    check.names = FALSE
  )
  res_f <- prep_uka_paired(cs_full, spec_cutoff = 0, control = "Control")
  expect_equal(res_f$fscore, 2.0)
  expect_equal(res_f$cell_line, "Treated")
})

test_that("prep_uka_paired flips sign when control is on the left of the contrast", {
  uka <- data.frame(
    `x.contrast` = c("Control vs Treated", "Treated2 vs Control"),
    `x.Kinase Name` = c("K1", "K2"),
    `x.Median Kinase Statistic` = c(1.0, 1.0),
    `x.Mean Specificity Score` = c(2.0, 2.0),
    check.names = FALSE
  )
  res <- prep_uka_paired(uka, spec_cutoff = 0, control = "Control", cs = FALSE)

  expect_equal(res$LogFC[res$cell_line == "Treated"], -1.0) # control on left -> sign flipped
  expect_equal(res$LogFC[res$cell_line == "Treated2"], 1.0) # control on right -> unchanged
})

test_that("prep_uka_paired keeps only comparisons that involve the control", {
  uka <- data.frame(
    `x.contrast` = c("DrugA vs DMSO", "DrugB vs DrugC", "DMSO vs DrugD"),
    `x.Kinase Name` = c("K1", "K2", "K3"),
    `x.Median Kinase Statistic` = c(1.0, 1.0, 1.0),
    `x.Mean Specificity Score` = c(2.0, 2.0, 2.0),
    check.names = FALSE
  )
  res <- prep_uka_paired(uka, spec_cutoff = 0, control = "DMSO", cs = FALSE)
  expect_setequal(res$cell_line, c("DrugA", "DrugD")) # "DrugB vs DrugC" dropped
})

test_that("prep_sens computes fold-change vs control and can restrict to best_drug_per_target", {
  sens <- data.frame(
    CELL_LINE_NAME = c("Control", "Control", "A", "A"),
    TARGET_1 = c("T1", "T2", "T1", "T2"),
    LN_IC50 = c(1.0, 2.0, 3.0, 1.5),
    DRUG_ID = c("D1", "D2", "D1", "D2")
  )
  res <- prep_sens(sens, control = "Control")
  expect_equal(res$LogFC[res$uniprotname == "T1"], 2.0) # 3.0 - 1.0
  expect_equal(res$LogFC[res$uniprotname == "T2"], -0.5) # 1.5 - 2.0
  expect_false("Control" %in% res$cell_line)

  res_restricted <- prep_sens(sens, control = "Control", best_drug_per_target = data.frame(DRUG_ID = "D1"))
  expect_equal(unique(res_restricted$uniprotname), "T1")
})
