test_that("make_golden_score dispatches to kinase-only when sens is NULL, full when it's given", {
  called <- NULL
  local_mocked_bindings(
    make_golden_score_kinase = function(uka, ...) called <<- "kinase",
    make_golden_score_full = function(uka, sens, ...) called <<- "full",
    .package = "networkScore"
  )

  expect_output(make_golden_score(uka = "x"), "kinase-only")
  expect_equal(called, "kinase")

  expect_output(make_golden_score(uka = "x", sens = "y"), "full")
  expect_equal(called, "full")
})
