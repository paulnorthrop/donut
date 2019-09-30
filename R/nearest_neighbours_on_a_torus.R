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
#' @return An object (a list) of class \code{c("nnt", "donut")} containing the
#'   following components.
#'   \item{nn.idx}{An \eqn{N} by \eqn{d} integer matrix of the \code{k}
#'     nearest neighbour indices, i.e. the rows of \code{data}.}
#'   \item{nn.dists}{An \eqn{N} by \eqn{d} numeric matrix of the \code{k}
#'     nearest neighbour distances.}
#'   \item{data, query, k}{The input arguments \code{data}, \code{query} and
#'     \code{k}.}
#'   \item{fn, torus, ranges, method}{The input arguments \code{fn},
#'     \code{torus}, \code{ranges} and \code{method}.}
#'   \item{call}{The call to \code{spm}.}
#' @seealso \code{\link[RANN:nn2]{RANN::nn2}},
#'   \code{\link[RANN.L1:nn2]{RANN.L1::nn2}},
#'   \code{\link[nabor:knn]{nabor::knn}}: nearest neigbour searches.
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
    res <- c(res, list(data = data, query = query, k = k, call = Call))
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
  # Check that all data and query are in the appropriate ranges
  # Check that ranges are in order (or make them in order?)
  # ................................
  #
  if (!(method %in% 1:2)) {
    stop("method must be equal to 1 or 2")
  }
  if (method == 1) {
    res <- method1_function(data, query, k, torus, ranges, fn, ...)
  } else {
    res <- method2_function(data, query, k, torus, ranges, fn, ...)
  }
  res <- c(res, list(data = data, query = query, k = k, fn = fn,
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
#'   when \code{ncol(x$data)} = 1 or 2.
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
#'   argument \code{col}.
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
