#' Maximum number of a trous levels for a signal of length n
#'
#' The B3 kernel at level j spans 2^(j+2) - 3 points. With
#' \code{criterion = "half"} the last level must not span more than half the
#' contour (a level that wraps around the whole closed contour carries no
#' local information); \code{"full"} allows up to the whole contour.
#'
#' @param n Integer. Signal length.
#' @param criterion \code{"half"} (default) or \code{"full"}.
#' @return Integer. For n = 512: 6 with \code{"half"}, 7 with \code{"full"}.
#' @export
max_levels <- function(n, criterion = c("half", "full")) {
  criterion <- match.arg(criterion)
  lim <- if (criterion == "half") n / 2 else n
  as.integer(floor(log2(lim + 3)) - 2L)
}

#' A trous wavelet decomposition with the B3-spline kernel
#'
#' Non-decimated (stationary) decomposition of a periodic signal. Each level
#' smooths the previous approximation with the B3 kernel (1, 4, 6, 4, 1)/16
#' dilated by 2^(j-1) (the "holes"), with periodic boundary; the detail is
#' the difference between consecutive approximations. Shift invariant and
#' zero phase, which is why it suits closed contours.
#'
#' @param x Numeric vector: the centroid-distance signature \code{r},
#'   periodic. A coordinate matrix is rejected.
#' @param J Integer. Number of levels. \code{NULL} (default) applies the
#'   half-turn criterion, \code{\link{max_levels}}; explicit values above the
#'   full-contour bound are rejected.
#'
#' @return A list with \code{c}, the \code{J + 1} approximations
#'   (\code{c[[1]]} is \code{x}), and \code{w}, the \code{J} details,
#'   \code{w[[j]] = c[[j]] - c[[j + 1]]}, so that
#'   \code{x = c[[J + 1]] + sum(w)} exactly.
#' @export
atrous_b3 <- function(x, J = NULL) {

  if (!is.numeric(x) || !is.null(dim(x)))
    stop("atrous_b3(): 'x' must be a numeric vector (the signature r), ",
         "not a matrix.", call. = FALSE)
  if (anyNA(x))
    stop("atrous_b3(): 'x' contains NA values.", call. = FALSE)
  x <- as.numeric(x)                 # drops attributes such as 'centroid'
  n <- length(x)

  if (is.null(J)) {
    J <- max_levels(n)               # half-turn criterion
  } else {
    if (length(J) != 1L || !is.numeric(J) || is.na(J) || J != round(J) || J < 1L)
      stop("atrous_b3(): 'J' must be a single positive integer or NULL.",
           call. = FALSE)
    J <- as.integer(J)
    if (J > max_levels(n, "full"))
      stop(sprintf(paste0("atrous_b3(): J = %d is too deep for n = %d: level %d ",
                          "spans %d points. Maximum %d (full contour), %d ",
                          "(half-turn)."),
                   J, n, J, 2^(J + 2) - 3, max_levels(n, "full"), max_levels(n)),
           call. = FALSE)
  }

  h <- c(1, 4, 6, 4, 1) / 16
  i <- seq_len(n) - 1L

  smooth_step <- function(sig, step) {       # periodic convolution, dilated by step
    y <- numeric(n)
    for (k in seq_along(h))
      y <- y + h[k] * sig[(i + (k - 3L) * step) %% n + 1L]
    y
  }

  c_list <- vector("list", J + 1L)
  c_list[[1L]] <- x
  for (j in seq_len(J))
    c_list[[j + 1L]] <- smooth_step(c_list[[j]], 2^(j - 1L))
  w_list <- lapply(seq_len(J), function(j) c_list[[j]] - c_list[[j + 1L]])

  list(c = c_list, w = w_list)
}
