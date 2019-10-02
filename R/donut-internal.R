#' Internal donut functions
#'
#' Internal donut functions
#' @details
#' These functions are not intended to be called by the user.
#' @name donut-internal
#' @keywords internal
NULL

# ============================ Function for method 1 ======================== #

#' @keywords internal
#' @rdname donut-internal
method1_function <- function(data, query, k, torus, ranges, fn, ...) {
  # The number of variables to be wrapped
  nt <- length(torus)
  # The respective ranges of each of these variables
  diffs <- apply(ranges, 1, diff)
  # 3^nt factorial design
  x <- fac3(nt)
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
  # Check for repeated indices in the rows of nnres$nn.idx
  dups <- apply(nnres$nn.idx, 1, anyDuplicated)
  # If any row has duplicated indices then call method2_function() for
  # the rows in query for which this is true
  if (any(dups > 0)) {
    which_dup <- which(dups > 0)
    res <- method2_function(data, query[which_dup, , drop = FALSE], k, torus,
                            ranges, fn, ...)
    nnres$nn.idx[which_dup, ] <- res$nn.idx
    nnres$nn.dists[which_dup, ] <- res$nn.dists
  }
  return(nnres)
}

# ============================ Function for method 2 ======================== #

#' @keywords internal
#' @rdname donut-internal
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

# ============================ 3^d factorial design ======================== #

#' @keywords internal
#' @rdname donut-internal
fac3 <- function(d) {
  # Creates a design matrix for a 3^d full factorial design
  #
  # Args:
  #   d : a positive integer scalar
  #
  # Returns:
  #   a 3^d by d integer matrix containing -1, 0, 1 for the levels of the d
  #   factors.  The first column varies fastest, the second column second
  #   fastest etc.
  rep_fn <- function(x) {
    rep(-1L:1L, times = 3L ^ x, each = 3L ^ (d - 1L - x))
  }
  return(vapply((d - 1L):0, rep_fn, numeric(3L ^ d)))
}
