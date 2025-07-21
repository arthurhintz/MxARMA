testthat::test_that("pdf_maxwell works correctly", {
  expect_snapshot(dmax(1, c(2, 5)))
  expect_snapshot(dmax(2, 1))
  expect_snapshot(dmax(5, c(1.5, 0.9, 0.8)))
  expect_snapshot(dmax(10, c(1, 0.912)))
  expect_snapshot(dmax(20, c(3.6789)))
  expect_snapshot(dmax(50, c(10, 20)))
})
