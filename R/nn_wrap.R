# ================================== nnt ==================================== #

#' Nearest Neighbour Search
#'
#' Uses a user-supplied function to find the \code{k} nearest neighbours of
#' specified points in a dataset, adding the option to wrap certain variables
#' on a torus.
#'
#' @param data An \eqn{M} by \eqn{d} numeric matrix or data frame.  Each of the
#'   \eqn{M} rows contains a \eqn{d}-dimensional observation.
#' @param query An \eqn{N} by \eqn{d} numeric matrix or data frame.  Each row
#'   contains an \eqn{d}-dimensional point that will be queried against
#'   \code{data}.
#' @param k An integer scalar.  The number of nearest neighbours, of the
#'   points in the rows of \code{query}, to find.
#' @param fn The function with which to calculate the nearest neighbours.
#'   The syntax of this function must be \code{fn(data, query, k, ...)}.
#'   The default is \code{RANN::nn2}.  Other possibilities are
#'   \code{RANN.L1:nn2} and \code{nabor::knn}.
#' @param torus An integer vector with element in
#'   \{1, ..., \code{ncol(data)}\}.
#' @param ranges A \code{length(torus)} by \code{2} numeric matrix.
#'   Row \code{i} gives the range of variation of the variable indexed by
#'   \code{torus[i]}. \code{ranges[i, 1]} and \code{ranges[i, 2]}
#'   are equivalent values of the variable, such as 0 degrees and 360 degrees.
#'   If \code{length(torus)} = 1 then \code{ranges} may be a vector of length
#'   2.
#' @param method An integer scalar, equal to 1 or 2.  See \strong{Details}.
#' @param ... Further arguments to be passed to \code{fn}.
#' @details
#'   If \code{method = 1} then the data are partially replicated, arranged
#'   around the original data in a way that wraps the variables in \code{torus} on their respective
#'   ranges in \code{ranges}.  Then \code{fn} is called using this replicated
#'   dataset as the argument \code{data}.  If \code{method = 2} then the
#'   following approach is used for the point in each row in \code{query}.
#'   The data indexed by \code{torus} are shifted (and wrapped) so that the
#'   point is located at the respective midpoints of \code{ranges}.
#'   This is only be an efficient approach if the number of points in
#'   \code{query} are is small.
#'
#'   If \code{torus} is missing then \code{fn} is called using
#'   \code{fn(data = data, query = query, k = k, ...)}, so that a call to
#'   \code{nnt} is equivalent to a call to the function chosen by \code{fn}.
#' @return An object (a list) of class \code{"nnt"} containing the following
#'   components.
#'   \item{nn.idx}{An \eqn{N} by \eqn{d} integer matrix of the \code{k}
#'     nearest neighbour indices, i.e. the rows of \code{data}.}
#'   \item{nn.dists}{An \eqn{N} by \eqn{d} numeric matrix of the \code{k}
#'     nearest neighbour distances.}
#'   \item{data, query}{The input arguments \code{data} and \code{query}.}
#'   \item{call}{The call to \code{spm}.}
#' @seealso \code{\link[RANN:nn2]{RANN::nn2}},
#'   \code{\link[RANN.L1:nn2]{RANN.L1::nn2}},
#'   \code{\link[nabor:knn]{nabor::knn}}: nearest neigbour searchess.
#' @examples
#' set.seed(20092019)
#' x1 <- c(0.2, 0.9)
#' x2 <- c(0.4, 0.7)
#' x <- cbind(x1, x2)
#' edge <- x
#' res <- nnt(x, edge)
#' ranges <- c(0, 1)
#' res <- nnt(x, edge, torus = 1, ranges = ranges)
#'
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
nnt <- function(data, query = data, k = min(10, nrow(data)),
                fn = RANN::nn2, torus, ranges, method = 1, ...) {
  Call <- match.call(expand.dots = TRUE)
  # Do the search and add data and query to the returned object
  if (missing(torus)) {
    res <- fn(data = data, query = query, k = k, ...)
    res <- c(res, list(data = data, query = query, call = Call))
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
  if (method %in% 1:2) {
    stop("method must be equal to 1 or 2")
  }
  if (method == 1) {
    res <- method1_function(data, query, k, torus, ranges, fn, ...)
  } else {
    res <- method2_function(data, query, k, torus, ranges, fn, ...)
  }
  res <- c(res, list(data = data, query = query, call = Call))
  class(res) <- c("nnt")
  return(res)
}

# ============================ Function for method 1 ======================== #

# Remove repeated indices
method1_function <- function(data, query, k, torus, ranges, fn, ...) {
  # The number of variables to be wrapped
  nt <- length(torus)
  # The respective ranges of each of these variables
  diffs <- apply(ranges, 1, diff)
  # 3^nt factorial design.  Need nVars = 1 to deal with the nt = 1 case
  x <- as.matrix(AlgDesign::gen.factorial(rep(3, nt), nVars = 1))
  # Multiply the -1s, 0s and 1s by the relevant values of diffs to surround the
  # original data by replicates in all directions
  x <- sweep(x, 2, diffs, "*")
  # Midpoints of each range
  mids <- rowMeans(ranges)
  # Function to replicate the data around the original data
  # Only use the variables to be wrapped
  # Each i corresponds to a different replicate
  tdata <- data[, torus, drop = FALSE]
  myfn <- function(i) {
    # Figure out which rows of tdata should be included
    # If x[i, j] < 0 then tdata[i, j] must be > mids[j]
    # If x[i, j] > 0 then tdata[i, j] must be < mids[j]
    # If x[i, j] = 0 then tdata[i, j] is unconstrained
    sgn <- sign(x[i, ])
    data_comp <- sweep(tdata, 2, sgn, "*")
    comp <- ifelse(sgn == 0, Inf, sgn * mids)
    cond <- sweep(data_comp, 2, comp, "<")
    which_rows <- apply(cond, 1, all)
    # Create the replicated data, add the other variables (if any)
    subdata <- tdata[which_rows, , drop = FALSE]
    subdata <- sweep(subdata, 2, x[i, ], "+")
    restdata <- data[which_rows, -torus, drop = FALSE]
    # Add the row numbers, so that we can return the original indices later
    idx <- which(which_rows)
    return(list(subdata = subdata, restdata = restdata, idx = idx))
  }
  # Return a list containing each replicated dataset (and the original)
  # and the indices of original data for each observation in rep_data
  # rep_data is an nrow(data)*(1+3^(nt-1)) by nt by matrix
  res <- vapply(1:nrow(x), myfn, list(subdata = 0, idx = 0, restdata = 0))
  rep_data <- do.call(rbind, res[1, ])
  rest_data <- do.call(rbind, res[2, ])
  # Replicate data
  big_data <- matrix(NA, nrow(rep_data), ncol(data), byrow = TRUE)
  big_data[, torus] <- rep_data
  big_data[, -torus] <- rest_data
  # Do the search using rep_data and the original query values
  nnres <- fn(data = big_data, query = query, k = k, ...)
  # Return to the indices that relate to the orginal data
  idx <- do.call(c, res[3, ])
  nnres$nn.idx <- apply(nnres$nn.idx, 1:2, function(x) idx[x])
  return(nnres)
}

# ============================ Function for method 2 ======================== #

method2_function <- function(data, query, k, torus, ranges, fn, ...) {
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
    nnres <- fn(data = ndata, query = nquery, k = k, ...)
    return(nnres)
  }
  # Call by_query() for each
  res <- vapply(1:nrow(query), by_query, list(nn.idx = 0, nn.dists = 0))
  nn.idx <- do.call(rbind, res[1, ])
  nn.dists <- do.call(rbind, res[2, ])
  res <- list(nn.idx = nn.idx, nn.dists = nn.dists)
  return(res)
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
