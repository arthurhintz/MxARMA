testthat::test_that("quantile function works correctly", {
  expect_snapshot(qmax(0.344, 2))
  expect_snapshot(qmax(0.78, 1))
  expect_snapshot(qmax(0.99999, 0.9))
  expect_snapshot(qmax(0.001, 0.112))
  expect_snapshot(qmax(0.294, 3.6789))
  expect_snapshot(qmax(0.5, 20))
})
