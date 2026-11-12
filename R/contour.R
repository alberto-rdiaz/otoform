#' Resample a closed contour at equal arc-length intervals
#'
#' Redistributes the vertices of a closed contour so that consecutive points
#' are equally spaced along the curve. The starting position can be shifted
#' with \code{s0}, which allows the sampling grid to be anchored at an
#' arbitrary point.
#'
#' @param xy A two-column matrix of contour coordinates.
#' @param N Number of output points.
#' @param s0 Arc-length offset of the first sample. Defaults to zero.
#'
#' @return An N-by-2 matrix, with the total perimeter stored as the
#'   attribute \code{perimeter}.
#' @export
#'
resample_arc <- function(xy, N, s0 = 0) {
  stopifnot(is.matrix(xy), ncol(xy) == 2L, nrow(xy) >= 3L,
            is.numeric(N), length(N) == 1L, N >= 3L)

  if (isTRUE(all.equal(xy[1L, ], xy[nrow(xy), ], check.attributes = FALSE))) {
    xy <- xy[-nrow(xy), , drop = FALSE]
  }

  xyc <- rbind(xy, xy[1L, ])

  seg <- sqrt(rowSums(diff(xyc)^2))
  if (any(seg == 0)) {
    xyc <- xyc[c(TRUE, seg > 0), , drop = FALSE]
    if (nrow(xyc) < 4L) stop("contorno degenerado tras filtrar segmentos nulos")
    seg <- sqrt(rowSums(diff(xyc)^2))
  }

  s <- c(0, cumsum(seg))
  P <- s[length(s)]

  t <- (s0 + seq(0, P, length.out = N + 1L)[-(N + 1L)]) %% P

  out <- cbind(x = stats::approx(s, xyc[, 1L], xout = t)$y,
               y = stats::approx(s, xyc[, 2L], xout = t)$y)

  attr(out, "perimeter") <- P
  out
}

#' Smooth a closed contour with a binomial kernel
#'
#' Applies a circular \[1/4, 1/2, 1/4\] kernel k times. This is equivalent to
#' a single convolution with a binomial kernel of order 2k, with standard
#' deviation sqrt(k/2) in units of vertex spacing.
#'
#' @param M A two-column matrix of contour coordinates.
#' @param k Number of smoothing passes. Zero leaves the contour unchanged.
#'
#' @return A matrix of the same dimensions as \code{M}.
#' @export
#'
smoothout <- function(M, k) {
  n  <- nrow(M)
  i  <- seq_len(n)
  im <- c(n, i[-n])          # indice del punto anterior (circular)
  ip <- c(i[-1L], 1L)        # indice del punto siguiente (circular)

  for (j in seq_len(k)) {
    M <- 0.5 * M + 0.25 * M[im, ] + 0.25 * M[ip, ]
  }
  M
}

#' Centre a contour and rotate it onto its principal axis
#'
#' Translates the contour so that the area centroid of the enclosed region is
#' at the origin and rotates it so that the major principal axis of the region
#' lies on the x axis. Centre and axis come from the polygon's area moments
#' (\code{\link{polygon_moments}}), so neither depends on how densely the
#' boundary is sampled, and the centre is the same point used later by
#' \code{\link{centroid_distance}}.
#'
#' @details
#' A principal axis is a line, not a direction: \code{theta} and
#' \code{theta + pi} describe the same axis. The angle is restricted to
#' \code{(-pi/2, pi/2]}, so the contour is never turned by more than a quarter
#' turn and whatever pointed right in the image still points right afterwards.
#' The 180-degree ambiguity is therefore resolved by the photographic protocol
#' of the package: every otolith is photographed with the rostrum to the
#' right. The function does not check this.
#'
#' If the two principal moments are (nearly) equal the axis is undefined; the
#' contour is returned centred but unrotated, with \code{status = "isotropic"}
#' and \code{theta = NA}.
#'
#' @param xy Numeric matrix, N x 2, vertices of a closed contour in pixels.
#' @param tol Numeric. Anisotropy below which the contour is treated as
#'   isotropic. Anisotropy is \code{(l1 - l2) / (l1 + l2)} for the principal
#'   moments \code{l1 >= l2}: 0 for a circle, 1 for a segment.
#'
#' @return The centred and rotated N x 2 matrix, with attributes
#'   \code{theta} (angle removed, radians; \code{NA} if isotropic),
#'   \code{center} (area centroid subtracted, input units),
#'   \code{anisotropy} and \code{status} (\code{"ok"} or \code{"isotropic"}).
#'
#' @export
rotation <- function(xy, tol = 1e-3) {

  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L)
    stop("rotation(): 'xy' must be a numeric matrix with 2 columns.",
         call. = FALSE)
  if (nrow(xy) < 3L)
    stop("rotation(): 'xy' needs at least 3 vertices.", call. = FALSE)
  if (anyNA(xy))
    stop("rotation(): 'xy' contains NA values.", call. = FALSE)

  pm     <- polygon_moments(xy)            # stops on zero area
  center <- pm$centroid
  xc     <- sweep(xy, 2L, center)          # keeps dimnames
  sxx <- pm$mxx; syy <- pm$myy; sxy <- pm$mxy

  mid <- (sxx + syy) / 2
  amp <- sqrt(((sxx - syy) / 2)^2 + sxy^2)
  anisotropy <- amp / mid                  # (l1 - l2) / (l1 + l2)

  if (anisotropy < tol) {
    out    <- xc
    theta  <- NA_real_
    status <- "isotropic"
  } else {
    theta <- 0.5 * atan2(2 * sxy, sxx - syy)   # in (-pi/2, pi/2]
    R <- matrix(c(cos(theta), sin(theta),
                  -sin(theta), cos(theta)), 2L, 2L, byrow = TRUE)
    out <- xc %*% t(R)
    dimnames(out) <- dimnames(xy)
    status <- "ok"
  }

  attr(out, "theta")      <- theta
  attr(out, "center")     <- center
  attr(out, "anisotropy") <- anisotropy
  attr(out, "status")     <- status
  out
}

#' Resample a contour with the rightmost point anchored at index 1
#'
#' Locates the rightmost point of the contour to sub-pixel accuracy by
#' fitting a parabola to the x coordinate near the discrete maximum, then
#' resamples the contour so that the first point falls exactly there.
#'
#' Assumes the otolith is oriented with the rostrum pointing right.
#'
#' @param xy A two-column matrix of contour coordinates.
#' @param N Number of output points. Defaults to 512.
#'
#' @return An N-by-2 matrix with attributes \code{status}, \code{s_star} and
#'   \code{delta}. Status is "ok" when anchoring succeeded.
#' @export
#'
anchor_rostrum <- function(xy, N = 512L) {
  g <- resample_arc(xy, N, s0 = 0)
  P <- attr(g, "perimeter")
  h <- P / N

  k  <- which.max(g[, 1L])
  km <- if (k == 1L) N else k - 1L
  kp <- if (k == N) 1L else k + 1L

  x0 <- g[km, 1L]; x1 <- g[k, 1L]; x2 <- g[kp, 1L]
  den <- x0 - 2 * x1 + x2

  if (den >= 0) return(structure(g, status = "pico_no_concavo"))

  delta <- (x0 - x2) / (2 * den)
  if (abs(delta) > 0.5) return(structure(g, status = "delta_fuera_rango"))

  s_star <- ((k - 1L + delta) * h) %% P

  structure(resample_arc(xy, N, s0 = s_star),
            status = "ok", s_star = s_star, delta = delta)
}

#' Normalise a closed contour to unit area
#'
#' Rescales a closed contour so that the area it encloses equals 1. Scaling by
#' \code{sqrt(A)} instead of by the maximum centroid distance makes the
#' normalisation depend on the whole shape rather than on a single vertex
#' (typically the rostrum tip), the vertex most exposed to extraction noise.
#'
#' @param xy Numeric matrix, N x 2, vertices of a closed contour in pixel
#'   units, typically the output of \code{\link{anchor_rostrum}}. Attributes
#'   are carried over: \code{status} and \code{delta} unchanged,
#'   \code{perimeter} and \code{s_star} rescaled to the new units.
#'
#' @return The rescaled N x 2 matrix, with one additional attribute,
#'   \code{area}: the enclosed area of the input in squared pixels, so the
#'   original scale is not lost.
#'
#' @details This is the last geometric step of the pipeline, after smoothing
#' and anchoring: every earlier step is scale-equivariant, and
#' \code{smoothout} slightly shrinks the enclosed area, so normalising earlier
#' would leave a shape-dependent residual in the scale.
#'
#' @seealso \code{\link{signed_area}}
#' @export
normalize_area <- function(xy) {

  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L)
    stop("normalize_area(): 'xy' must be a numeric matrix with 2 columns.",
         call. = FALSE)
  if (nrow(xy) < 3L)
    stop("normalize_area(): 'xy' needs at least 3 vertices.", call. = FALSE)
  if (anyNA(xy))
    stop("normalize_area(): 'xy' contains NA values.", call. = FALSE)

  A <- abs(signed_area(xy))
  if (!is.finite(A) || A <= 0)
    stop(sprintf("normalize_area(): enclosed area is %g; degenerate contour.",
                 A), call. = FALSE)

  s   <- sqrt(A)
  out <- xy / s   # arithmetic keeps dim, dimnames and all other attributes

  # Attributes in pixel units must change units with the coordinates.
  if (!is.null(attr(out, "perimeter")))
    attr(out, "perimeter") <- attr(xy, "perimeter") / s
  if (!is.null(attr(out, "s_star")))
    attr(out, "s_star") <- attr(xy, "s_star") / s
  attr(out, "area") <- A

  out
}

#' Area, centroid and second central moments of a closed polygon
#'
#' Shoelace-type sums over the edges. Every quantity refers to the enclosed
#' region, not to the vertices, so none depends on how densely the boundary
#' is sampled. The result is the same for clockwise and counter-clockwise
#' traversal. Assumes validated input.
#'
#' @param xy Numeric matrix, N x 2, vertices of a simple closed polygon.
#'
#' @return A list: \code{area} (positive), \code{centroid} (\code{c(x, y)}),
#'   and the central second moments \code{mxx} (\eqn{\int (x - C_x)^2 dA}),
#'   \code{myy} (\eqn{\int (y - C_y)^2 dA}) and \code{mxy}
#'   (\eqn{\int (x - C_x)(y - C_y) dA}).
#' @keywords internal
polygon_moments <- function(xy) {

  x  <- xy[, 1L]
  y  <- xy[, 2L]
  xn <- c(x[-1L], x[1L])
  yn <- c(y[-1L], y[1L])
  cross <- x * yn - xn * y

  A <- signed_area(xy)
  if (!is.finite(A) || A == 0)
    stop("polygon_moments(): zero or undefined area.", call. = FALSE)
  s <- sign(A)          # orientation-independent moments
  A <- abs(A)

  cx <- s * sum((x + xn) * cross) / (6 * A)
  cy <- s * sum((y + yn) * cross) / (6 * A)

  mxx <- s * sum((x * x + x * xn + xn * xn) * cross) / 12 - A * cx * cx
  myy <- s * sum((y * y + y * yn + yn * yn) * cross) / 12 - A * cy * cy
  mxy <- s * sum((x * yn + 2 * x * y + 2 * xn * yn + xn * y) * cross) / 24 -
    A * cx * cy

  list(area = A, centroid = c(cx, cy), mxx = mxx, myy = myy, mxy = mxy)
}

#' Centroid of a closed polygon
#'
#' The area centroid, from \code{\link{polygon_moments}}. Assumes validated
#' input.
#' @param xy Numeric matrix, N x 2, vertices of a closed polygon.
#' @return Numeric vector \code{c(x, y)}.
#' @keywords internal
polygon_centroid <- function(xy) polygon_moments(xy)$centroid

#' Centroid-distance signature of a closed contour
#'
#' Distance from the centroid of the enclosed region to each vertex, in the
#' units of \code{xy}. This is the 1-D signal the wavelet decomposition works
#' on. The centroid is that of the polygon's area (see
#' \code{\link{polygon_centroid}}), not the mean of the vertices: on real
#' otolith contours the two differ by tens of pixels because the rostrum
#' contributes much boundary for little area.
#'
#' @param xy Numeric matrix, N x 2, vertices of a closed contour, typically
#'   the output of \code{\link{normalize_area}}. No rescaling is applied here.
#'
#' @return Numeric vector \code{r} of length N, with attribute
#'   \code{centroid}: the centroid \code{c(x, y)} used.
#'
#' @seealso \code{\link{normalize_area}}, \code{\link{atrous_b3}}
#' @export
centroid_distance <- function(xy) {

  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L)
    stop("centroid_distance(): 'xy' must be a numeric matrix with 2 columns.",
         call. = FALSE)
  if (nrow(xy) < 3L)
    stop("centroid_distance(): 'xy' needs at least 3 vertices.", call. = FALSE)
  if (anyNA(xy))
    stop("centroid_distance(): 'xy' contains NA values.", call. = FALSE)
  if (signed_area(xy) == 0)
    stop("centroid_distance(): degenerate contour (zero area).", call. = FALSE)

  cen <- polygon_centroid(xy)
  r   <- sqrt((xy[, 1L] - cen[1L])^2 + (xy[, 2L] - cen[2L])^2)
  attr(r, "centroid") <- cen
  r
}

#' Rotate a contour so that its rightmost point comes first
#'
#' Cyclically shifts a closed contour so that the point with the maximum x
#' coordinate becomes the first element. Assumes the otolith is oriented
#' with the rostrum pointing right.
#'
#' @param out A list with elements X and Y, the contour coordinates.
#'
#' @return A list with elements x and y, cyclically shifted.
#' @export
#'
coordsSort <- function(out) {
  k <- which.max(out[[1]])
  n <- length(out[[1]])
  i <- c(k:n, seq_len(k - 1L))
  list(x = out[[1]][i],
       y = out[[2]][i])
}
