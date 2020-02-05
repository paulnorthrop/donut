#' donut: Nearest Neighbour Search with Variables on a Torus
#'
#' Finds the \code{k} nearest neighbours in a dataset of specified points,
#' adding the option to wrap certain variables on a torus.  The user chooses
#' the algorithm to use to find the nearest neighbours.
#'
#' @details The function \code{\link{nnt}} performs the nearest neighbour
#'   search.  There is also a rudimentary plot method: \code{\link{plot.nnt}}.
#'
#'   The default algorithm is that provided by the function
#'   \code{\link[RANN]{nn2}} in the \code{\link[RANN]{RANN-package}}.
#'   Another possibility is the \code{\link[nabor]{knn}} function in the
#'   \code{\link[nabor]{nabor-package}}.
#'
#'   See \code{vignette("donut-vignette", package = "donut")} for an
#'   overview of the package.
#' @references Arya, S., Mount, D., Kemp, S. E. and Jefferis, G. (2019)
#'   RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using L2
#'   Metric. R package version 2.6.1.
#'   \url{https://CRAN.R-project.org/package=RANN}
#' @references Elseberg J., Magnenat S., Siegwart R., Nuchter, A. (2012)
#'   Comparison of nearest-neighbor-search strategies and implementations for
#'   efficient shape registration. \emph{Journal of Software Engineering for
#'   Robotics (JOSER)}, \strong{3}(1), 2-12
#'   \url{https://CRAN.R-project.org/package=nabor}
#' @seealso \code{\link{nnt}} for nearest neighbour with some variables
#'   wrapped on a torus.
#' @seealso \code{\link{plot.nnt}} plot method for objects returned from
#'   \code{\link{nnt}} (1 and 2 dimensional data only).
#' @docType package
#' @name donut
NULL
