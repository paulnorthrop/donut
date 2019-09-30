context("Mop up")

# Check that RANN is available
got_RANN <- requireNamespace("RANN", quietly = TRUE)

if (got_RANN) {
  set.seed(29092019)
  x1 <- runif(100, 0, 2 * pi)
  x2 <- runif(100, 0, 3)
  x3 <- runif(100, -1, 1)

  # Call nnt() with torus but missing ranges
  x <- cbind(x1, x2, x3)
  test_that("nnt: no ranges", {
    the_error <- "ranges must be supplied"
    testthat::expect_error(nnt(x, x, torus = 1, method = 1), the_error)
  })
}
