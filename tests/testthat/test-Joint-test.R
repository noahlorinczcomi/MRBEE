test_that("Joint.test returns a data frame with Chi2 and P columns", {
  set.seed(42)
  n <- 50
  p <- 3
  bZ <- matrix(rnorm(n * p), n, p)
  RZ <- diag(p)
  res <- suppressMessages(capture.output(
    out <- Joint.test(bZ, RZ)
  ))
  expect_s3_class(out, "data.frame")
  expect_named(out, c("Chi2", "P"))
  expect_equal(nrow(out), n)
})

test_that("Joint.test Chi2 values are non-negative and P-values are in [0, 1]", {
  set.seed(7)
  n <- 100
  p <- 2
  bZ <- matrix(rnorm(n * p), n, p)
  RZ <- diag(p)
  capture.output(out <- Joint.test(bZ, RZ))
  expect_true(all(out$Chi2 >= 0))
  expect_true(all(out$P >= 0 & out$P <= 1))
})

test_that("Joint.test returns large Chi2 for strong signals", {
  set.seed(1)
  n <- 10
  p <- 2
  bZ <- matrix(10, n, p)
  RZ <- diag(p)
  capture.output(out <- Joint.test(bZ, RZ))
  expect_true(all(out$P < 0.001))
})
