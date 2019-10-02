# ================================== nnt ==================================== #

#' Nearest Neighbour Search with Variables on a Torus
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
#' @param torus An integer vector with element(s) in
#'   \{1, ..., \code{ncol(data)}\}.  The corresponding variables are wrapped
#'   on the corresponding range gives in \code{ranges}.
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
#'   dataset as the argument \code{data}.  If \code{k} is large and/or
#'   \code{data} is a sparse dataset then it is possible that a single
#'   observation contributes more than once to a set of nearest neighbours,
#'   which is incorrect. If this occurs then \code{nnt} uses method 2 to
#'   correct the offending rows in \code{nn.idx} and \code{nn.dists} in the
#'   returned list object.
#'
#'   If \code{method = 2} then the
#'   following approach is used for the point in each row in \code{query}.
#'   The data indexed by \code{torus} are shifted (and wrapped) so that the
#'   point is located at the respective midpoints of \code{ranges}.
#'   Method 2 is efficient only if the number of points in \code{query} is
#'   small.
#'
#'   If \code{torus} is missing then \code{fn} is called using
#'   \code{fn(data = data, query = query, k = k, ...)}, so that a call to
#'   \code{nnt} is equivalent to a call to the function chosen by \code{fn}.
#' @return An object (a list) of class \code{c("nnt", "donut")} containing the
#'   following components.
#'   \item{nn.idx}{An \eqn{N} by \eqn{d} integer matrix of the \code{k}
#'     nearest neighbour indices, i.e. the rows of \code{data}.}
#'   \item{nn.dists}{An \eqn{N} by \eqn{d} numeric matrix of the \code{k}
#'     nearest neighbour distances.}
#'   \item{data, query, k, fn}{The arguments \code{data}, \code{query},
#'     \code{k} and \code{fn} (in fact \code{substitute(fn)}).}
#'   \item{torus, ranges, method}{If \code{torus} is supplied, the
#'     arguments \code{torus}, \code{ranges} and \code{method}.}
#'   \item{call}{The call to \code{spm}.}
#' @seealso \code{\link[RANN:nn2]{RANN::nn2}},
#'   \code{\link[RANN.L1:nn2]{RANN.L1::nn2}},
#'   \code{\link[nabor:knn]{nabor::knn}}: nearest neighbour searches.
#' @references Arya, S., Mount, D., Kemp, S. E. and Jefferis, G. (2019)
#'   RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using L2
#'   Metric. R package version 2.6.1.
#'   \url{https://CRAN.R-project.org/package=RANN}
#' @references Arya, S., Mount, D., Kemp, S. E., Jefferis, G. and Muller,
#'   K. (2018)  RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using
#'   L1 Metric. R package version 2.5.2.
#'   \url{https://CRAN.R-project.org/package=RANN.L1}
#' @references Elseberg J., Magnenat S., Siegwart R., Nuchter, A. (2012)
#'   Comparison of nearest-neighbor-search strategies and implementations for
#'   efficient shape registration. \emph{Journal of Software Engineering for
#'   Robotics (JOSER)}, \strong{3}(1), 2-12
#'   \url{https://CRAN.R-project.org/package=nabor}
#' @seealso \code{\link{plot.nnt}} plot method for objects returned from
#'   \code{\link{nnt}} (1 and 2 dimensional data only).
#' @examples
#' got_RANN <- requireNamespace("RANN", quietly = TRUE)
#' got_RANN.L1 <- requireNamespace("RANN.L1", quietly = TRUE)
#' got_nabor <- requireNamespace("nabor", quietly = TRUE)
#'
#' set.seed(20092019)
#' # 2D example from the RANN:nn2 documentation (L2 metric)
#' x1 <- runif(100, 0, 2 * pi)
#' x2 <- runif(100, 0, 3)
#' DATA <- data.frame(x1, x2)
#' if (got_RANN) {
#'   nearest <- nnt(DATA, DATA)
#' }
#'
#' # Suppose that x1 should be wrapped
#' ranges1 <- c(0, 2 * pi)
#' query1 <- rbind(c(6, 1.3), c(2 * pi, 3), c(3, 1.5), c(4, 0))
#' if (got_RANN) {
#'   res1 <- nnt(DATA, query1, k = 8, torus = 1, ranges = ranges1)
#'   plot(res1, ylim = c(0, 3))
#' }
#'
#' # Suppose that x1 and x2 should be wrapped
#' ranges2 <- rbind(c(0, 2 * pi), c(0, 3))
#' query2 <- rbind(c(6, 1.3), c(2 * pi, 3), c(3, 1.5), c(4, 0))
#' if (got_RANN) {
#'   res2 <- nnt(DATA, query2, k = 8, torus = 1:2, ranges = ranges2)
#'   plot(res2)
#' }
#'
#' # Use nabor::knn (L2 metric) instead of RANN::nn2
#' if (got_nabor) {
#'   res3 <- nnt(DATA, query2, k = 8, fn = nabor::knn, torus = 1:2,
#'               ranges = ranges2)
#'   plot(res3)
#' }
#'
#' # Use RANN.L1::nn2 (L1 metric)
#' if (got_RANN.L1) {
#'   res4 <- nnt(DATA, query2, k = 8, fn = RANN.L1::nn2, torus = 1:2,
#'               ranges = ranges2)
#'   plot(res4)
#' }
#'
#' # 1D example
#' ranges <- c(0, 2 * pi)
#' query <- c(4, 0.1)
#' if (got_RANN) {
#'   res <- nnt(x1, query, torus = 1, ranges = ranges, method = 1)
#'   plot(res)
#' }
#' @export
nnt <- function(data, query = data, k = min(10, nrow(data)),
                fn = RANN::nn2, torus, ranges, method = 1, ...) {
  Call <- match.call(expand.dots = TRUE)
  # Make data and query matrices
  data <- as.matrix(data)
  query <- as.matrix(query)
  # Do the search and add data and query to the returned object
  if (missing(torus)) {
    res <- fn(data = data, query = query, k = k, ...)
    res <- c(res, list(data = data, query = query, k = k, fn = substitute(fn),
                       call = Call))
    class(res) <- c("nnt", "donut")
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
  # Sort the values in the rows of ranges, from smallest to largest
  ranges <- t(apply(ranges, 1, sort))
  # Check that all data and query are in the appropriate ranges
  check_ranges <- function(i, x) {
    return(any(x[, torus[i]] < ranges[i, 1] | x[, torus[i]] > ranges[i, 2]))
  }
  check_data <- vapply(1:length(torus), check_ranges, TRUE, x = data)
  check_query <- vapply(1:length(torus), check_ranges, TRUE, x = query)
  if (any(check_data)) {
    stop("value(s) in 'data' are outside the corresponding range in 'ranges'")
  }
  if (any(check_query)) {
    stop("value(s) in 'query' are outside the corresponding range in 'ranges'")
  }
  # Check method
  if (!(method %in% 1:2)) {
    stop("method must be equal to 1 or 2")
  }
  if (method == 1) {
    res <- method1_function(data, query, k, torus, ranges, fn, ...)
  } else {
    res <- method2_function(data, query, k, torus, ranges, fn, ...)
  }
  res <- c(res, list(data = data, query = query, k = k, fn = substitute(fn),
                     torus = torus, ranges = ranges, method = method,
                     call = Call))
  class(res) <- c("nnt", "donut")
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
#'   or \code{\link[graphics]{points}}.
#' @details This function is only applicable in 1 or 2 dimensions, that is,
#'   when \code{ncol(x$data)} = 1 or 2.  It provides a visual check that the
#'   wrapping of variables is working as intended, in cases where the
#'   number of query points, that is, \code{nrow(x$query)} is small
#'   enough that sets of nearest neighbours do not overlap much.
#'
#'   If \code{ncol(x$data)} = 1 then the index of each observation is plotted
#'   against its value, using a plotting character \code{pch = 1}.  A vertical
#'   line is superimposed at each value in \code{x$query} and the \code{x$k$}
#'   nearest neighbours of each line are colour-coded.
#'
#'   If \code{ncol(x$data)} = 2 then \code{x$data[, 2]} is plotted against
#'   \code{x$data[, 1]}, using a plotting character \code{pch = 1}.  Each point
#'   in \code{x$query} is plotted with a cross and the \code{x$k$}
#'   nearest neighbours of each point are colour-coded.
#'
#'   Colours of the lines/crosses and nearest neighbour points can be set sing an
#'   argument \code{col}.  If a variable is wrapped then the default plotting
#'   limits are set using the corresponding values in \code{x$ranges}.
#' @return Nothing is returned.
#' @seealso \code{\link{nnt}} for nearest neighbour with some variables
#'   wrapped on a torus.
#' @section Examples:
#' See the examples in \code{\link{nnt}}.
#' @export
plot.nnt <- function(x, ...) {
  if (!inherits(x, "donut")) {
    stop("use only with \"donut\" objects")
  }
  ncov <- ncol(x$data)
  if (ncov > 2) {
    stop("The plot method works for up to 2 covariates only")
  }
  if (is.null(colnames(x$data))) {
    colnames(x$data) <- paste0("X", 1:ncol(x$data))
  }
  my_plot <- function(x, ..., pch, lwd, col, xlim = my_xlim, ylim = my_ylim) {
    graphics::plot(x, ..., lwd = 1, col = "black", xlim = xlim, ylim = ylim)
  }
  my_points <- function(x, ..., pch = 16, col = "red", lwd) {
    graphics::points(x, ..., pch = pch, col = col, lwd = 1)
  }
  user_args <- list(...)
  nquery <- nrow(x$query)
  if (is.null(user_args$col)) {
    user_args$col <- 1 + 1:nquery
  }
  my_xlim <- my_ylim <- NULL
  if (ncov == 2) {
    if (!is.null(x$torus)) {
      if (1 %in% x$torus) {
        my_xlim <- x$ranges[1, ]
      }
      if (2 %in% x$torus) {
        my_ylim <- x$ranges[2, ]
      }
    }
    my_plot(x$data, ...)
    for (i in 1:nquery) {
      i_user_args <- user_args
      i_user_args$col <- user_args$col[i]
      for_my_points <- c(list(x = x$data[x$nn.idx[i, ], , drop = FALSE]),
                         i_user_args)
      do.call(my_points, for_my_points)
    }
    for_points <- list(x = x$query, col = user_args$col, pch = "x")
    do.call(graphics::points, for_points)
  } else {
    if (!is.null(x$torus)) {
      my_xlim <- x$ranges
    }
    plot_data <- cbind(x$data, index = 1:nrow(x$data))
    my_plot(plot_data, ...)
    graphics::abline(v = x$query, col = user_args$col)
    for (i in 1:nquery) {
      i_user_args <- user_args
      i_user_args$col <- user_args$col[i]
      plot_data <- cbind(x$data[x$nn.idx[i, ], , drop = FALSE], x$nn.idx[i, ])
      for_my_points <- c(list(x = plot_data), i_user_args)
      do.call(my_points, for_my_points)
    }
  }
  return(invisible())
}
