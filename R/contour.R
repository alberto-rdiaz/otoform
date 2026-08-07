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

  out <- cbind(x = approx(s, xyc[, 1L], xout = t)$y,
               y = approx(s, xyc[, 2L], xout = t)$y)

  attr(out, "perimeter") <- P
  out
}
