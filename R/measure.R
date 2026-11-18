# R/measure.R

#' Shape measurements of an aligned contour
#'
#' Classical size and shape descriptors of a closed contour, computed from
#' the coordinates (never from attributes). On the output of
#' \code{\link{anchor_rostrum}} the major axis lies on x, so \code{length}
#' and \code{width} are the extents along the principal axes.
#'
#' @param xy Numeric matrix, N x 2, the aligned contour \code{points} in
#'   pixels.
#' @param cal Pixels per millimetre. \code{NULL} (default) keeps pixel units;
#'   a single positive number converts lengths (/ cal) and areas (/ cal^2)
#'   to millimetres. Dimensionless descriptors do not change.
#'
#' @return A named numeric vector: \code{area}, \code{perimeter},
#'   \code{length}, \code{width}, \code{roundness} (\eqn{4 \pi A / P^2}, 1 for
#'   a circle), \code{aspect} (\code{length / width}) and
#'   \code{rectangularity} (\code{area / (length * width)}), with attribute
#'   \code{unit} (\code{"px"} or \code{"mm"}).
#' @export
measure_contour <- function(xy, cal = NULL) {

  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L)
    stop("measure_contour(): 'xy' must be a numeric matrix with 2 columns.",
         call. = FALSE)
  if (nrow(xy) < 3L)
    stop("measure_contour(): 'xy' needs at least 3 vertices.", call. = FALSE)
  if (anyNA(xy))
    stop("measure_contour(): 'xy' contains NA values.", call. = FALSE)
  if (!is.null(cal) && (length(cal) != 1L || !is.numeric(cal) ||
                        is.na(cal) || cal <= 0))
    stop("measure_contour(): 'cal' must be NULL or a single positive number ",
         "(pixels per mm).", call. = FALSE)

  A  <- polygon_moments(xy)$area
  xn <- c(xy[-1L, 1L], xy[1L, 1L])
  yn <- c(xy[-1L, 2L], xy[1L, 2L])
  P  <- sum(sqrt((xn - xy[, 1L])^2 + (yn - xy[, 2L])^2))
  L  <- diff(range(xy[, 1L]))
  W  <- diff(range(xy[, 2L]))

  s <- if (is.null(cal)) 1 else cal
  out <- c(area           = A / s^2,
           perimeter      = P / s,
           length         = L / s,
           width          = W / s,
           roundness      = 4 * pi * A / P^2,
           aspect         = L / W,
           rectangularity = A / (L * W))
  attr(out, "unit") <- if (is.null(cal)) "px" else "mm"
  out
}
