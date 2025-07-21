testthat::test_that("cdf_maxwell works correctly", {
  expect_snapshot(pmax(1, c(2, 5)))
  expect_snapshot(pmax(2, c(1)))
  expect_snapshot(pmax(5, c(1.5, 0.9, 0.8)))
  expect_snapshot(pmax(10, c(1, 0.912)))
  expect_snapshot(pmax(20, c(3.6789)))
  expect_snapshot(pmax(50, c(10, 20)))
})
