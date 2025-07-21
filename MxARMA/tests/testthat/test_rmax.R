testthat::test_that("generating random variables setted", {
  set.seed(8)
  expect_snapshot(rmax(c(2, 5)))
  expect_snapshot(rmax(1))
  expect_snapshot(rmax(c(1.5, 0.9, 0.8)))
  expect_snapshot(rmax(c(1, 2, 3, 4)))
  expect_snapshot(rmax(100))
  expect_snapshot(rmax(c(10, 20)))
})
