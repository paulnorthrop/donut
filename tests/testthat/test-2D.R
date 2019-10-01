context("2D check by hand")

# Check that RANN is available
got_RANN <- requireNamespace("RANN", quietly = TRUE)

# Consider cases where method 1 (on its own) gives the wrong nearest neighbours
# because an observation, or pair of observations, contribute more than once.
# Check that an additional call using method 2 inside nnt() corrects this.

if (got_RANN) {
  x1 <- c(5, 5)
  x2 <- c(9, 3)
  x3 <- c(8, 9)
  x4 <- c(0, 2)
  x5 <- c(1, 10)
  x <- rbind(x1, x2, x3, x4, x5)
#  plot(x, xlim = c(0, 10), ylim = c(0, 10))

  dfn <- function(x, y) sqrt(sum((x - y) ^ 2))

  # Only wrap on variable 1
  n1d <- n1n <- matrix(NA, 5, 5)
  temp <- apply(cbind(x1, x2, x3, x4, x5), 2, dfn, y = x1)
  n1d[1, ] <- sort(temp)
  n1n[1, ] <- order(temp)
  x4d <- c(10, 2)
  x5d <- c(11, 10)
  temp <- apply(cbind(x1, x2, x3, x4d, x5d), 2, dfn, y = x2)
  n1d[2, ] <- sort(temp)
  n1n[2, ] <- order(temp)
  temp <- apply(cbind(x1, x2, x3, x4d, x5d), 2, dfn, y = x3)
  n1d[3, ] <- sort(temp)
  n1n[3, ] <- order(temp)
  x2d <- c(-1, 3)
  x3d <- c(-2, 9)
  temp <- apply(cbind(x1, x2d, x3d, x4, x5), 2, dfn, y = x4)
  n1d[4, ] <- sort(temp)
  n1n[4, ] <- order(temp)
  temp <- apply(cbind(x1, x2d, x3d, x4, x5), 2, dfn, y = x5)
  n1d[5, ] <- sort(temp)
  n1n[5, ] <- order(temp)
  res1 <- nnt(x, x, torus = 1, ranges = c(0, 10), method = 1)
  res2 <- nnt(x, x, torus = 1, ranges = c(0, 10), method = 2)
  test_that("Wrap on variable 1: indices vs method 1", {
    testthat::expect_equal(res1$nn.idx, n1n)
  })
  test_that("Wrap on variable 1: indices vs method 2", {
    testthat::expect_equal(res2$nn.idx, n1n)
  })
  test_that("Wrap on variable 1: distances vs method 1", {
    testthat::expect_equal(res1$nn.dists, n1d)
  })
  test_that("Wrap on variable 1: distances vs method 2", {
    testthat::expect_equal(res2$nn.dists, n1d)
  })

  # Only wrap on variable 2
  n2d <- n2n <- matrix(NA, 5, 5)
  temp <- apply(cbind(x1, x2, x3, x4, x5), 2, dfn, y = x1)
  n2d[1, ] <- sort(temp)
  n2n[1, ] <- order(temp)
  x3d <- c(8, -1)
  x5d <- c(1, 0)
  temp <- apply(cbind(x1, x2, x3d, x4, x5d), 2, dfn, y = x2)
  n2d[2, ] <- sort(temp)
  n2n[2, ] <- order(temp)
  x2d <- c(9, 13)
  x4d <- c(0, 12)
  temp <- apply(cbind(x1, x2d, x3, x4d, x5), 2, dfn, y = x3)
  n2d[3, ] <- sort(temp)
  n2n[3, ] <- order(temp)
  x3d <- c(8, -1)
  x5d <- c(1, 0)
  temp <- apply(cbind(x1, x2, x3d, x4, x5d), 2, dfn, y = x4)
  n2d[4, ] <- sort(temp)
  n2n[4, ] <- order(temp)
  temp <- apply(cbind(x1, x2d, x3, x4d, x5), 2, dfn, y = x5)
  n2d[5, ] <- sort(temp)
  n2n[5, ] <- order(temp)
  res1 <- nnt(x, x, torus = 2, ranges = c(0, 10), method = 1)
  res2 <- nnt(x, x, torus = 2, ranges = c(0, 10), method = 2)
  test_that("Wrap on variable 2: indices vs method 1", {
    testthat::expect_equal(res1$nn.idx, n2n)
  })
  test_that("Wrap on variable 2: indices vs method 2", {
    testthat::expect_equal(res2$nn.idx, n2n)
  })
  test_that("Wrap on variable 2: distances vs method 1", {
    testthat::expect_equal(res1$nn.dists, n2d)
  })
  test_that("Wrap on variable 2: distances vs method 2", {
    testthat::expect_equal(res2$nn.dists, n2d)
  })

  # Wrap on variables 1 and 2
  n2d <- n2n <- matrix(NA, 5, 5)
  temp <- apply(cbind(x1, x2, x3, x4, x5), 2, dfn, y = x1)
  n2d[1, ] <- sort(temp)
  n2n[1, ] <- order(temp)
  x3d <- c(8, -1)
  x4d <- c(10, 2)
  x5d <- c(11, 0)
  temp <- apply(cbind(x1, x2, x3d, x4d, x5d), 2, dfn, y = x2)
  n2d[2, ] <- sort(temp)
  n2n[2, ] <- order(temp)
  x2d <- c(9, 13)
  x4d <- c(10, 12)
  x5d <- c(11, 10)
  temp <- apply(cbind(x1, x2d, x3, x4d, x5d), 2, dfn, y = x3)
  n2d[3, ] <- sort(temp)
  n2n[3, ] <- order(temp)
  x2d <- c(-1, 3)
  x3d <- c(-2, -1)
  x5d <- c(1, 0)
  temp <- apply(cbind(x1, x2d, x3d, x4, x5d), 2, dfn, y = x4)
  n2d[4, ] <- sort(temp)
  n2n[4, ] <- order(temp)
  x2d <- c(-1, 13)
  x3d <- c(-2, 9)
  x4d <- c(0, 12)
  temp <- apply(cbind(x1, x2d, x3d, x4d, x5), 2, dfn, y = x5)
  n2d[5, ] <- sort(temp)
  n2n[5, ] <- order(temp)
  ranges <- rbind(c(0, 10), c(0, 10))
  res1 <- nnt(x, x, torus = 1:2, ranges = ranges, method = 1)
  res2 <- nnt(x, x, torus = 1:2, ranges = ranges, method = 2)
  test_that("Wrap on variable 2: indices vs method 1", {
    testthat::expect_equal(res1$nn.idx, n2n)
  })
  test_that("Wrap on variable 2: indices vs method 2", {
    testthat::expect_equal(res2$nn.idx, n2n)
  })
  test_that("Wrap on variable 2: distances vs method 1", {
    testthat::expect_equal(res1$nn.dists, n2d)
  })
  test_that("Wrap on variable 2: distances vs method 2", {
    testthat::expect_equal(res2$nn.dists, n2d)
  })
}

