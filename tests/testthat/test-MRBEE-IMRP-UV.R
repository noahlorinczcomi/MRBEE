set.seed(42)
n <- 100
bx <- rnorm(n, sd = 0.5)
bxse <- rep(0.1, n)
theta_true <- 0.5
by <- bx * theta_true + rnorm(n, sd = 0.1)
byse <- rep(0.1, n)
Rxy <- diag(2)

test_that("MRBEE.IMRP.UV returns a correctly structured list", {
  res <- MRBEE.IMRP.UV(by = by, bx = bx, byse = byse, bxse = bxse, Rxy = Rxy)
  expect_type(res, "list")
  expect_true(all(c("theta", "vartheta", "delta") %in% names(res)))
  expect_length(res$theta, 1)
  expect_true(is.numeric(res$vartheta))
  expect_true(res$vartheta > 0)
  expect_length(res$delta, n)
})

test_that("MRBEE.IMRP.UV bootstrap variance option returns sampling.theta", {
  res <- MRBEE.IMRP.UV(
    by = by, bx = bx, byse = byse, bxse = bxse, Rxy = Rxy,
    var.method = "bootstrap", sampling.time = 50
  )
  expect_true("sampling.theta" %in% names(res))
  expect_true(is.numeric(res$sampling.theta))
})
