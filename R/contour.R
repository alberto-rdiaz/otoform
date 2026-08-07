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
