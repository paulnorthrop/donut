context("Plot method")

# Check that RANN is available
got_RANN <- requireNamespace("RANN", quietly = TRUE)

if (got_RANN) {
  set.seed(29092019)
  x1 <- runif(100, 0, 2 * pi)
  x2 <- runif(100, 0, 3)
  x3 <- runif(100, -1, 1)

  # 3D
  x <- cbind(x1, x2, x3)
  res <- nnt(x, x, torus = 1, ranges = c(0, 2 * pi), method = 1)
  test_that("RANN: plot 3D", {
    the_error <- "The plot method works for up to 2 covariates only"
    testthat::expect_error(plot(res), the_error)
  })

  # 2D
  x <- cbind(x1, x2)
  ranges <- matrix(c(0, 0, 2 * pi, 3), 2, 2)
  query <- rbind(c(6, 0.1), c(3, 3))
  res <- nnt(x, query, torus = 1:2, ranges = ranges, method = 1)
  test_that("RANN: plot 2D", {
    testthat::expect_equal(plot(res), NULL)
  })

  # 1D
  ranges <- c(0, 2 * pi)
  query <- c(4, 0.1)
  res <- nnt(x1, query, torus = 1, ranges = ranges, method = 1)
  test_that("RANN: plot 1D", {
    testthat::expect_equal(plot(res), NULL)
  })
}
