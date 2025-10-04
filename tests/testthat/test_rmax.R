testthat::test_that("generating random variables setted", {
  set.seed(8)
  expect_snapshot(rmax(1, 0.5))
  expect_snapshot(rmax(2, 1))
  expect_snapshot(rmax(3, 1.5))
  expect_snapshot(rmax(4, 2))
})
