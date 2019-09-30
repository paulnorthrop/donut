context("No torus")

set.seed(20092019)
# Example from the RANN:nn2 documentation
x1 <- runif(100, 0, 2 * pi)
x2 <- runif(100, 0,3)
DATA <- data.frame(x1, x2)

# Check that RANN is available
got_RANN <- requireNamespace("RANN", quietly = TRUE)
# No wrapping
if (got_RANN) {
  # Use RANN::nn2 directly (this is the default fn)
  nearest <- RANN::nn2(DATA,DATA)
  # Use donut::nnt
  res <- nnt(DATA, DATA, fn = RANN::nn2)
  test_that("RANN: no wrapping, nn.idx", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
  test_that("RANN: no wrapping, nn.dists", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
}

# Check that RANN.L1 is available
got_RANN.L1 <- requireNamespace("RANN.L1", quietly = TRUE)
# No wrapping
if (got_RANN.L1) {
  # Use RANN.L1::nn2 directly
  nearest <- RANN.L1::nn2(DATA,DATA)
  # Use donut::nnt
  res <- nnt(DATA, DATA, fn = RANN.L1::nn2)
  test_that("RANN.L1: no wrapping, nn.idx", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
  test_that("RANN.L1: no wrapping, nn.dists", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
}

# Check that nabor is available
got_nabor <- requireNamespace("nabor", quietly = TRUE)
# No wrapping
if (got_nabor) {
  # Use nabor::knn directly
  nearest <- nabor::knn(DATA,DATA, k = 10)
  # Use donut::nnt
  res <- nnt(DATA, DATA, fn = nabor::knn)
  test_that("nabor: no wrapping, nn.idx", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
  test_that("nabor: no wrapping, nn.dists", {
    testthat::expect_equal(nearest$nn.idx, res$nn.idx)
  })
}
