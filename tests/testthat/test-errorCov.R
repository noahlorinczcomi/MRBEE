test_that("errorCov returns a symmetric matrix with ones on the diagonal", {
  set.seed(42)
  n <- 300
  p <- 3
  ZMatrix <- matrix(rnorm(n * p), n, p)
  res <- capture.output(
    out <- errorCov(ZMatrix, Zscore.cutoff = 0.5, subsampling.ratio = 0.5, subsampling.time = 20)
  )
  expect_equal(dim(out), c(p, p))
  expect_equal(out, t(out), tolerance = 1e-10)
  expect_equal(diag(out), rep(1, p), tolerance = 1e-10)
})

test_that("errorCov off-diagonal values are in (-1, 1)", {
  set.seed(99)
  n <- 300
  p <- 2
  ZMatrix <- matrix(rnorm(n * p), n, p)
  capture.output(out <- errorCov(ZMatrix, Zscore.cutoff = 0.5, subsampling.ratio = 0.5, subsampling.time = 20))
  expect_true(all(abs(out[lower.tri(out)]) < 1))
})
