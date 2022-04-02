#context("Repeated indices are deleted")

# Check that RANN is available
got_RANN <- requireNamespace("RANN", quietly = TRUE)

# Consider cases where method 1 (on its own) gives the wrong nearest neighbours
# because an observation, or pair of observations, contribute more than once.
# Check that an additional call using method 2 inside nnt() corrects this.

if (got_RANN) {
  x1 <- c(0.9, 0.9)
  x2 <- c(0.1, 0.9)
  x3 <- c(0.5, 0.1)
  x <- rbind(x1, x2, x3)

  # Only wrap on variable 1
  res1 <- nnt(x, x, torus = 1, ranges = c(0, 1), method = 1)
  res2 <- nnt(x, x, torus = 1, ranges = c(0, 1), method = 2)
  # call and method are different, but otherwise the results should be the same
  res1$call <- res2$call <- NULL
  res1$method <- res2$method <- NULL
  test_that("RANN: repeats 1D", {
    testthat::expect_equal(res1, res2)
  })

  # Wraps on variables 1 and 2
  x4 <- c(0.05, 0.05)
  x5 <- c(0.85, 0.15)
  x <- rbind(x1, x2, x3, x4, x5)
  ranges <- matrix(c(0, 0, 1, 1), 2, 2)
  query <- rbind(x, c(0.1, 0.5), c(1, 0.4))
  res1 <- nnt(x, query, torus = 1:2, ranges = ranges, method = 1)
  res2 <- nnt(x, query, torus = 1:2, ranges = ranges, method = 2)
  # call and method are different, but otherwise the results should be the same
  res1$call <- res2$call <- NULL
  res1$method <- res2$method <- NULL
  test_that("RANN: repeats 2D", {
    testthat::expect_equal(res1, res2)
  })
}

