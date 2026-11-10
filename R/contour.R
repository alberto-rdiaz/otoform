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

#' Centre and rotate a contour onto its principal axis
#'
#' Translates the contour to the origin and rotates it so that its principal
#' axis aligns with the x axis. The 180-degree ambiguity of the principal
#' axis is resolved by assuming the rostrum points right in the source image.
#'
#' @param M A two-column matrix of contour coordinates.
#'
#' @return A two-column matrix, centred and rotated, with the rotation angle
#'   stored as the attribute \code{theta}.
#' @export
#'
rotation <- function(M) {
  stopifnot(is.matrix(M), ncol(M) == 2L, nrow(M) >= 3L)

  # Centrar en el origen
  Mc <- cbind(M[, 1L] - mean(M[, 1L]),
              M[, 2L] - mean(M[, 2L]))

  # Eje principal: autovector asociado al autovalor mayor
  v <- eigen(stats::cov(Mc))$vectors[, 1L]
  theta <- atan2(v[2L], v[1L])

  R <- matrix(c( cos(theta), sin(theta),
                 -sin(theta), cos(theta)), 2L, 2L, byrow = TRUE)

  Mr <- Mc %*% t(R)

  # Resolver la ambiguedad de 180 grados:
  # el rostrum debe quedar en x positiva
  if (abs(min(Mr[, 1L])) > abs(max(Mr[, 1L]))) {
    Mr <- -Mr
    theta <- theta + pi
  }

  attr(Mr, "theta") <- theta
  Mr
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
