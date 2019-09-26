# ================================ nn2wrap ================================== #

#' Nearest Neighbour Search
#'
#' Uses \code{\link[RANN]{nn2}} to find the nearest neighbours in a dataset
#' specified points, adding the option to wrap certain variables on a torus.
#'
#' @param data An \eqn{M} by \eqn{d} numeric matrix or data frame.  Each of the
#'   \eqn{M} rows contains a \eqn{d}-dimensional observation.
#' @param query An \eqn{N} by \eqn{d} numeric matrix or data frame.  Each row
#'   contains an \eqn{d}-dimensional point that will be queried against
#'   \code{data}.
#' @param k An integer scalar.  The number of nearest neighbours, of the
#'   points in the rows of \code{query}, to find.
#' @param torus An integer vector with element in
#'   \{1, ..., \code{ncol(data)}\}.  If \code{torus} is missing then
#'   a call to \code{nnt} is equivalent to a call to the function chosen
#'   by \code{package}.
#' @param ranges A \code{length(torus)} by \code{2} numeric matrix.
#'   Row \code{i} gives the range of variation of the variable indexed by
#'   \code{torus[i]}. \code{ranges[i, 1]} and \code{ranges[i, 2]}
#'   are equivalent values of the variable, such as 0 degrees and 360 degrees.
#'   If \code{length(torus)} = 2 then \code{ranges} may be a vector of length
#'   2.
#' @param package The package to use calculate the nearest neighbours.
#'   One of \code{"RANN"}, \code{"RANN.L1"} or \code{"nabor"}.
#' @param method An integer scalar, equal to 1 or 2.
#' @param ... Further arguments to be passed to
#'   \code{\link[RANN:nn2]{RANN::nn2}},
#'   \code{\link[RANN.L1:nn2]{RANN.L1::nn2}} or
#'   \code{\link[nabor:knn]{nabor::knn}}
#' @details Add details
#' @return Return
#' @seealso \code{\link[RANN:nn2]{RANN::nn2}},
#'   \code{\link[RANN.L1:nn2]{RANN.L1::nn2}},
#'   \code{\link[nabor:knn]{nabor::knn}}: nearest neigbour searchess.
#' @examples
#' set.seed(20092019)
#' # Example from the RANN:nn2 documentation
#' x1 <- runif(100, 0, 2 * pi)
#' x2 <- runif(100, 0, 3)
#' DATA <- data.frame(x1, x2)
#' nearest <- nnt(DATA, DATA)
#'
#' # Suppose that x1 should be wrapped
#' ranges <- c(0, 2 * pi)
#' nearest <- nnt(DATA, DATA, torus = 1, ranges = ranges)
#' edge <- matrix(c(2 * pi, 1.5), 1, 2)
#' res <- nnt(DATA, edge, torus = 1, ranges = ranges)
#' plot(res, pch = 16)
#'
#' # Suppose that x1 and x2 should be wrapped
#' ranges <- rbind(c(0, 2 * pi), c(0, 3))
#' nearest <- nnt(DATA, DATA[1, ], torus = 1:2, ranges = ranges)
#' edge <- rbind(c(2 * pi, 1.5), c(2 * pi, 3))
#' res <- nnt(DATA, edge, torus = 1:2, ranges = ranges)
#' plot(res, pch = 16)
#'
#' n <- 100
#' pjn <- cbind(runif(n, 0, 2), runif(n, 0, 1))
#' ranges <- rbind(c(0, 2), c(0, 1))
#' edge <- rbind(c(1.0, 0.5), c(0.1, 0.1))
#' res <- nnt(pjn, edge, torus = 1:2, ranges = ranges, method = 2)
#' plot(res, pch = 16)
#'
#' y <- nshs[, "hs"]
#' x <- nshs[, c("season", "direction")]
#' query <- matrix(c(350, 0, 150, 360), 2, 2)
#' ranges <- matrix(c(0, 0, 360, 360), 2, 2)
#' res <- nnt(data = x, query = query, k = 100, torus = 1:2, ranges = ranges)
#' res <- nnt(data = x, query = as.matrix(x), k = 100, torus = 1:2, ranges = ranges)
#' plot(res)
#'
#' which_vals <- 1:5
#' y <- nshs[which_vals, "hs"]
#' x <- nshs[which_vals, c("season", "direction")]
#' query <- matrix(c(350, 0, 150, 360), 2, 2)
#' ranges <- matrix(c(0, 0, 360, 360), 2, 2)
#' res <- nnt(data = x, query = query, k = 2, torus = 1:2, ranges = ranges)
#' @export
nnt <- function(data, query = data, k = min(10, nrow(data)), torus, ranges,
                     package = c("RANN", "RANN.L1", "nabor"), method = 1,
                ...) {
  package <- match.arg(package)
  # Check that the chosen package is available
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package ", package, "is not installed", call. = FALSE)
  }
  # Set the function to be used
  which_fn <- switch(package,
                     RANN = RANN::nn2,
                     RANN.L1 = RANN.L1::nn2,
                     nabor = nabor::knn)
  # Do the search and add data and query to the returned object
  if (missing(torus)) {
    res <- which_fn(data = data, query = query, k = k, ...)
    res <- c(res, list(data = data, query = query))
    class(res) <- c("nnt")
    return(res)
  }
  if (missing(ranges)) {
    stop("ranges must be supplied")
  }
  if (length(torus) == 0L || mode(torus) != "numeric") {
    stop("'torus' must be a non-empty numeric vector")
  }
  if (!all(torus %in% 1:ncol(data))) {
    stop("All elements of 'torus' must be in 1:ncol(data)")
  }
  dim_ranges <- ncol(ranges)
  if (is.null(dim_ranges)) {
    dim(ranges) <- 1:2
  }
  if (!is.matrix(ranges) || ncol(ranges) != 2 ||
      nrow(ranges) != length(torus)) {
    stop("ranges not consistent with length(torus)")
  }
  # Make data and query matrices
  data <- as.matrix(data)
  query <- as.matrix(query)
  # Check that all data and query are in the appropriate ranges
  # Check that ranges are in order (or make them in order?)
  # ................................
  # We transform the data in data and query, columns in torus
  # Repeat for each row in query
  #
  if (method == 1) {
    res <- method1_function(data, query, k, torus, ranges, which_fn, ...)
  } else {
    res <- method2_function(data, query, k, torus, ranges, which_fn, ...)
  }
  res <- c(res, list(data = data, query = query))
  class(res) <- c("nnt")
  return(res)
}

# ============================ Function for method 1 ======================== #

method1_function <- function(data, query, k, torus, ranges, which_fn, ...) {
  # Midpoints, ranges, lower and upper limits variables to be wrapped
  mids <- rowMeans(ranges)
  diffs <- apply(ranges, 1, diff)
  low <- ranges[, 1]
  high <- ranges[, 2]
  # Function to find the nearest neighbours for the jth row of query, that is,
  # for the jth point of interest
  by_query <- function(j) {
    # Take a copy of the data and shift data for the variables to be wrapped
    # so that the point of interest is at the midpoint of the each variable
    ndata <- data
    shoof <- query[j, torus] - mids
    ndata[, torus] <- sweep(ndata[, torus, drop = FALSE], 2, shoof, "-")
    # The consequence of the shift is that data are moved below or above the
    # ranges of the variables. Wrap these data round to the parts of the
    # ranges that have been vacated data that have been shifted up (shoof is
    # negative) or down (shoof is positive).
    #
    # The amounts that we need to move these data are stored in diffs.
    # The direction that we move them depends on the sign of shoof (sgn):
    # we will add sgn * diffs to such data
    sgn <- sign(shoof)
    # Multiply data for that were shifted up by -1.  Then we can identify
    # values that lie above the range by whether they are less than high
    data_comp <- sweep(ndata[, torus, drop = FALSE], 2, sgn, "*")
    # If sgn =  1: data were shifted down and values < low are shifted up
    # If sgn = -1: data were shifted up and -values < -high are shifted down
    comp <- (sgn == 1) * low - (sgn == -1) * high
    # Find the values to move (for whom comp = TRUE) and then move (only) them
    cond <- sweep(data_comp, 2, comp, "<")
    ndata[, torus] <- ndata[, torus] + sweep(cond, 2, sgn * diffs, "*")
    # Move the point of interest to the midpoint
    nquery <- query[j, , drop = FALSE]
    nquery[, torus] <- mids
    # Do the search using the new data and the new query values
    nnres <- which_fn(data = ndata, query = nquery, k = k, ...)
    return(nnres)
  }
  # Call by_query() for each
  res <- vapply(1:nrow(query), by_query, list(nn.idx = 0, nn.dists = 0))
  nn.idx <- do.call(rbind, res[1, ])
  nn.dists <- do.call(rbind, res[2, ])
  res <- list(nn.idx = nn.idx, nn.dists = nn.dists)
  return(res)
}

# ============================ Function for method 2 ======================== #

# Do only for variables in torus
# Add variables not in torus
# Remove repeated indices
method2_function <- function(data, query, k, torus, ranges, which_fn, ...) {
  # The number of variables to be wrapped
  nt <- length(torus)
  # The respective ranges of each of these variables
  diffs <- apply(ranges, 1, diff)
  # 3^nt factorial design
  x <- as.matrix(AlgDesign::gen.factorial(rep(3, nt)))
  # Multiply the -1s, 0s and 1s by the relevant values of diffs to surround the
  # original data by replicates in all directions
  x <- sweep(x, 2, diffs, "*")
  # Midpoints of each range
  mids <- rowMeans(ranges)
  # Function to replicate the data around the original data
  myfn <- function(i) {
    # Use only a subset of the data
    # If x[i, j] < 0 then data[i, j] must be > mids[j]
    # If x[i, j] > 0 then data[i, j] must be < mids[j]
    # If x[i, j] = 0 then data[i, j] is unconstrained
    sgn <- sign(x[i, ])
    data_comp <- sweep(data, 2, sgn, "*")
    comp <- ifelse(sgn == 0, Inf, sgn * mids)
    cond <- sweep(data_comp, 2, comp, "<")
    which_data <- apply(cond, 1, all)
    subdata <- data[which_data, , drop = FALSE]
    subdata <- sweep(subdata, 2, x[i, ], "+")
    idx <- which(which_data)
    return(list(subdata = subdata, idx = idx))
  }
  # Return a list containing each replicated dataset (and the original)
  # and the indices of original data for each observation in rep_data
  # rep_data is an nrow(data)*(1+3^(nt-1)) by nt by matrix
  res <- vapply(1:nrow(x), myfn, list(subdata = 0, idx = 0))
  rep_data <- do.call(rbind, res[1, ])
  idx <- do.call(c, res[2, ])
  # Do the search using rep_data and the original query values
  nnres <- which_fn(data = rep_data, query = query, k = k, ...)
  # Return to the indices that relate to the orginal data
  nnres$nn.idx <- apply(nnres$nn.idx, 1:2, function(x) idx[x])
  return(nnres)
}

# =========================== Plot nearest neighbours ======================= #

#' Plot diagnostics for an nnt object
#'
#' \code{plot} method for an objects of class \code{c("nnt")}.
#'
#' @param x an object of class \code{c("nnt")}, a result of
#'   a call to \code{\link{nnt}}.
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}},
#'    \code{\link[graphics]{lines}} or \code{\link[graphics]{points}}.
#' @return Nothing is returned.
#' @seealso \code{\link{set_val_data}}: to set validation data.
#' @section Examples:
#' See the examples in \code{\link{nnt}}.
#' @export
plot.nnt <- function(x, ...) {
  if (!inherits(x, "nnt")) {
    stop("use only with \"nnt\" objects")
  }
  ncov <- ncol(x$data)
  if (ncov > 2) {
    stop("The plot method works for up to 2 covariates only")
  }
  my_plot <- function(x, ..., lwd, col) {
    graphics::plot(x, ..., lwd = 1, col = "black")
  }
  my_points <- function(x, ..., col = "red", lwd) {
    graphics::points(x, ..., col = col, lwd = 1)
  }
  if (ncov == 2) {
    # Plot covariate positions: validation data in red
    my_plot(x$data, ...)
    for (i in 1:nrow(x$query)) {
      print(x$data[x$nn.idx[i, ], ])
      my_points(x$data[x$nn.idx[i, ], , drop = FALSE], ...)
      # Add the (circular) limits of the kernel
      theta <- seq(0, 2 * pi, len = 100)
      r <- max(x$nn.dists[i, ])
      xvals <- r * cos(theta) + x$query[i, 1]
      yvals <- r * sin(theta) + x$query[i, 2]
      graphics::lines(xvals, yvals, ...)
      graphics::points(x$query, pch = "x")
    }
  } else {
    # Plot covariate positions: validation data in red
    my_plot(cbind(x$x, x$y), ...)
    my_points(x$val_data[, c(2, 1)], ...)
    # Add the limits of the kernel
    v <- c(x$x0 - x$val_bw, x$x0, x$x0 + x$val_bw)
    graphics::abline(v = v, ...)
  }
  return(invisible())
}
