context("No torus")

# Check that RANN is available
got_RANN <- requireNamespace("RANN", quietly = TRUE)

# No wrapping

if (got_RANN) {
  set.seed(20092019)
  # Example from the RANN:nn2 documentation
  x1 <- runif(100, 0, 2*pi)
  x2 <- runif(100, 0,3)
  DATA <- data.frame(x1, x2)
  # Use RANN::nn2 directly
  nearest <- RANN::nn2(DATA,DATA)
  # Use donut::nnt
  res <- nnt(DATA, DATA)
  res$data <- NULL
  res$call <- NULL
  res$query <- NULL
  test_that("RANN: no wrapping, nn.idx", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
  test_that("RANN: no wrapping, nn.dists", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
}

