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

  # Call nnt() with non-numeric torus
  test_that("nnt: no ranges", {
    the_error <- "'torus' must be a non-empty numeric vector"
    testthat::expect_error(nnt(x, x, torus = "1", ranges = c(0, 1),
                               method = 1), the_error)
  })

  # Call nnt() with a torus that is not in {1, ..., ncol(data)}
  test_that("nnt: no ranges", {
    testthat::expect_error(nnt(x, x, torus = c(1:2, 5), ranges = c(0, 1),
                               method = 1))
  })

  # Call nnt() with a ranges that is not consistent with torus
  test_that("nnt: no ranges", {
    the_error <- "ranges not consistent with length(torus)"
    testthat::expect_error(nnt(x, x, torus = 1:2, ranges = c(0, 1),
                               method = 1))
  })
  # Call nnt() with an inappropriate value of method
  test_that("nnt: no ranges", {
    the_error <- "method must be equal to 1 or 2"
    testthat::expect_error(nnt(x, x, torus = 1, ranges = c(0, 1),
                               method = 3))
  })
}
