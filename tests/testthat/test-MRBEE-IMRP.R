set.seed(42)
n <- 100
p <- 2
bX <- matrix(rnorm(n * p, sd = 0.5), n, p)
bXse <- matrix(0.1, n, p)
theta_true <- c(0.3, -0.2)
by <- drop(bX %*% theta_true) + rnorm(n, sd = 0.1)
byse <- rep(0.1, n)
Rxy <- diag(p + 1)

test_that("MRBEE.IMRP returns a correctly structured list", {
  res <- MRBEE.IMRP(by = by, bX = bX, byse = byse, bXse = bXse, Rxy = Rxy)
  expect_type(res, "list")
  expect_named(res, c("theta", "covtheta", "delta"))
  expect_length(res$theta, p)
  expect_equal(dim(res$covtheta), c(p, p))
  expect_length(res$delta, n)
})

test_that("MRBEE.IMRP covariance matrix is symmetric positive semi-definite", {
  res <- MRBEE.IMRP(by = by, bX = bX, byse = byse, bXse = bXse, Rxy = Rxy)
  cov <- res$covtheta
  expect_equal(cov, t(cov), tolerance = 1e-10)
  expect_true(all(eigen(cov)$values >= -1e-10))
})

test_that("MRBEE.IMRP.Egger returns list with intercept in theta", {
  res <- MRBEE.IMRP.Egger(by = by, bX = bX, byse = byse, bXse = bXse, Rxy = Rxy)
  expect_type(res, "list")
  expect_named(res, c("theta", "covtheta", "delta"))
  expect_length(res$theta, p + 1)
  expect_equal(dim(res$covtheta), c(p + 1, p + 1))
  expect_length(res$delta, n)
})
