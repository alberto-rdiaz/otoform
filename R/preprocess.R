#' Convert an RGB image to luminance
#'
#' Applies the Rec.709 luminance coefficients (0.2126, 0.7152, 0.0722). A plain
#' channel average is not used: it over-weights the blue channel and flattens
#' the otolith/background contrast.
#'
#' @param img A 3D array (height x width x 3) with values in \code{[0, 1]},
#'   as returned by \code{jpeg::readJPEG}. Exactly three channels are
#'   required and they are interpreted as R, G and B in that order.
#'
#' @return A numeric matrix (height x width) with values in \code{[0, 1]}.
#'   The result is \strong{continuous}, not quantised.
#'
#' @details
#' The output is deliberately left unquantised. The Rec.709 weights are not
#' dyadic, so combining three 8-bit channels yields tens of thousands of
#' distinct luma values, not 256. Choosing a number of histogram levels
#' \code{L} is a separate methodological decision and belongs to
#' \code{\link{binarize}}: colour conversion and histogram binning are
#' justified on different grounds and should not be coupled.
#'
#' @seealso \code{\link{binarize}}, \code{\link{otsu}}
#'
#' @export
to_gray <- function(img) {

  if (!is.array(img) || !is.numeric(img) || length(dim(img)) != 3L)
    stop("to_gray(): 'img' must be a 3D array (height x width x 3). Greyscale ",
         "JPEGs are read as 2D matrices and carry no colour information.",
         call. = FALSE)

  if (dim(img)[3L] != 3L)
    stop(sprintf("to_gray(): 'img' has %d channel(s); RGB requires exactly 3.",
                 dim(img)[3L]), call. = FALSE)

  if (anyNA(img))
    stop("to_gray(): 'img' contains NA values in the RGB channels.",
         call. = FALSE)

  rng <- range(img)
  if (rng[1L] < 0 || rng[2L] > 1)
    stop(sprintf("to_gray(): values outside [0, 1]: observed [%g, %g]. ",
                 rng[1L], rng[2L]),
         "Some readers return [0, 255]; rescale before calling.",
         call. = FALSE)

  0.2126 * img[,,1L] + 0.7152 * img[,,2L] + 0.0722 * img[,,3L]
}

#' Binarise a luminance image
#'
#' Quantises the luma to \code{L} levels, applies a threshold (either chosen
#' automatically by Otsu's method or supplied by the user) and normalises the
#' polarity so that \code{TRUE} always marks the object.
#'
#' @param g Numeric matrix of luminance values in \code{[0, 1]}, as returned
#'   by \code{\link{to_gray}}.
#' @param method Either \code{"otsu"} (default) or \code{"fixed"}.
#' @param threshold Only with \code{method = "fixed"}: an integer in
#'   \code{[0, L - 2]}, expressed on the \strong{quantised scale}, not on
#'   \code{[0, 1]}. It carries the same meaning as the value returned by
#'   \code{\link{otsu}}: the last level belonging to the low class. Must be
#'   \code{NULL} when \code{method = "otsu"}.
#' @param L Integer. Number of histogram levels. Defaults to 256.
#'
#' @return A list with:
#'   \describe{
#'     \item{bw}{Logical matrix, \code{TRUE} marking the object. \code{NULL}
#'       for degenerate images.}
#'     \item{threshold}{The threshold actually applied, on the quantised
#'       scale. \code{NA_integer_} for degenerate images.}
#'     \item{separability}{Otsu's \code{eta*} evaluated at the applied
#'       threshold. Defined for both methods, which allows the loss incurred
#'       by a fixed threshold to be quantified against the optimum.}
#'     \item{inverted}{Logical. \code{TRUE} if the mask had to be flipped.}
#'     \item{method}{The method used.}
#'     \item{L}{Number of levels used.}
#'     \item{status}{\code{"ok"} or \code{"degenerate"}.}
#'   }
#'
#' @details
#' \strong{Scale.} Pixels are compared on the quantised scale
#' (\code{g_q > threshold}), never by rescaling the threshold back to
#' \code{[0, 1]}. The two are not equivalent: pixels whose \code{g * (L - 1)}
#' falls in \code{[k - 0.5, k)} round to \code{k} and would change class
#' depending on where the comparison is made. That half-level is a band of
#' pixels along the otolith edge, precisely where \code{Conte} traces.
#'
#' \strong{Comparison.} The test is \code{>} and not \code{>=}, because the
#' threshold is the last level of the low class.
#'
#' \strong{Polarity.} Otsu separates low from high levels but cannot know
#' which group is the object. Polarity is normalised using an invariant of the
#' photographic setup: the background touches the image border, a centred
#' otolith does not. If most of the frame comes out \code{TRUE} the mask is
#' flipped and \code{inverted} records it. A sudden rise in the proportion of
#' inverted frames means something changed in the setup.
#'
#' @seealso \code{\link{to_gray}}, \code{\link{otsu}}
#'
#' @export
binarize <- function(g, method = c("otsu", "fixed"), threshold = NULL,
                     L = 256L) {

  method <- match.arg(method)

  if (length(L) != 1L || !is.numeric(L) || is.na(L) || L != round(L) || L < 2L)
    stop("binarize(): 'L' must be a single integer >= 2.", call. = FALSE)
  L <- as.integer(L)

  if (!is.matrix(g) || !is.numeric(g))
    stop("binarize(): 'g' must be a numeric matrix.", call. = FALSE)

  if (anyNA(g))
    stop("binarize(): 'g' contains NA values.", call. = FALSE)

  rng <- range(g)
  if (rng[1L] < 0 || rng[2L] > 1)
    stop(sprintf("binarize(): 'g' outside [0, 1]: observed [%g, %g]. ",
                 rng[1L], rng[2L]),
         "Expected the output of to_gray().", call. = FALSE)

  # An argument silently ignored is a bug waiting to happen: it is an error.
  if (method == "otsu" && !is.null(threshold))
    stop("binarize(): 'threshold' does not apply with method = \"otsu\". ",
         "Use method = \"fixed\" to set it by hand.", call. = FALSE)

  if (method == "fixed") {
    if (is.null(threshold))
      stop("binarize(): method = \"fixed\" requires 'threshold'.",
           call. = FALSE)
    if (length(threshold) != 1L || !is.numeric(threshold) ||
        is.na(threshold) || threshold != round(threshold))
      stop("binarize(): 'threshold' must be a single integer.", call. = FALSE)
    if (threshold < 0L || threshold > L - 2L)
      stop(sprintf("binarize(): 'threshold' = %g outside [0, %d]. ",
                   threshold, L - 2L),
           "It is given on the quantised scale, not on [0, 1]: for a cut at ",
           sprintf("%.2f use %d.", 0.5, floor(0.5 * (L - 1L) - 0.5)),
           call. = FALSE)
  }

  # The only quantisation point in the whole pipeline.
  g_q <- round(g * (L - 1L))

  if (method == "otsu") {
    res <- otsu(g_q, L = L)
    thr <- res$threshold
    sep <- res$separability
    sts <- res$status
  } else {
    thr <- as.integer(threshold)
    st  <- otsu_stats(g_q, L)
    if (st$sigma_T == 0) {
      sep <- NA_real_
      sts <- "degenerate"
    } else {
      sep <- st$sigma_B[thr + 1L]
      sep <- if (is.finite(sep)) sep / st$sigma_T else NA_real_
      sts <- "ok"
    }
  }

  if (sts == "degenerate")
    return(list(bw = NULL, threshold = NA_integer_, separability = NA_real_,
                inverted = NA, method = method, L = L,
                status = "degenerate"))

  bw <- g_q > thr

  # Corners are counted twice; irrelevant for a mean over thousands of border
  # pixels.
  brd      <- c(bw[1L, ], bw[nrow(bw), ], bw[, 1L], bw[, ncol(bw)])
  inverted <- mean(brd) > 0.5
  if (inverted) bw <- !bw

  list(bw           = bw,
       threshold    = thr,
       separability = sep,
       inverted     = inverted,
       method       = method,
       L            = L,
       status       = "ok")
}

#' Compute the signed area of a closed polygon (shoelace formula)
#'
#' The sign of the result indicates the traversal direction: positive for
#' counter-clockwise, negative for clockwise. Used internally to enforce a
#' consistent contour orientation.
#'
#' @param M A 2D matrix of coordinates.
#'
#' @return A single numeric value. Positive if counter-clockwise.
#' @keywords internal
#'
signed_area <- function(M) {
  sum(M[,1L] * c(M[,2L][-1], M[,2L][1]) - c(M[,1L][-1], M[,1L][1]) * M[,2L]) / 2
}

#' Trace the contour of a binary object (Moore boundary tracing)
#'
#' Follows the outer boundary of a binary shape starting from a point found
#' by scanning inward from the image edge.
#'
#' @param I A 2D binary matrix (0s and 1s), typically the output of
#'   \code{binarize}.
#'
#' @return A two-column matrix of contour coordinates, in counter-clockwise
#'   order.
#' @export
#'
Conte <- function(I) {
  x <- (c(floor(dim(I)[1]/2), floor(dim(I)[2]/2)))
  x[1] <- dim(I)[1] - x[1]

  while (x[2] > 1 && abs(I[x[1], x[2]] - I[x[1], (x[2]-1)]) < 0.1) {
    x[2] <- x[2] - 1
  }

  a <- 1
  M <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),2,8,byrow=TRUE)
  M <- cbind(M[,8], M, M[,1])
  X <- 0; Y <- 0
  x1 <- x[1]; x2 <- x[2]
  SS <- NA; S <- 6

  while ((any(c(X[a], Y[a]) != c(x1, x2)) | length(X) < 3)) {
    if (abs(I[x[1]+M[1,S+1], x[2]+M[2,S+1]] - I[x[1],x[2]]) < 0.1) {
      a <- a+1; X[a] <- x[1]; Y[a] <- x[2]; x <- x+M[,S+1]
      SS[a] <- S+1; S <- (S+7) %% 8
    } else if (abs(I[x[1]+M[1,S+2], x[2]+M[2,S+2]] - I[x[1],x[2]]) < 0.1) {
      a <- a+1; X[a] <- x[1]; Y[a] <- x[2]; x <- x+M[,S+2]
      SS[a] <- S+2; S <- (S+7) %% 8
    } else if (abs(I[x[1]+M[1,(S+3)], x[2]+M[2,(S+3)]] - I[x[1],x[2]]) < 0.1) {
      a <- a+1; X[a] <- x[1]; Y[a] <- x[2]; x <- x+M[,(S+3)]
      SS[a] <- S+3; S <- (S+7) %% 8
    } else S <- (S+1) %% 8
  }

  out <- cbind(x = Y[-1], y = (dim(I)[1] - X)[-1])

  if (signed_area(out) < 0) {
    out <- out[rev(seq_len(nrow(out))), , drop = FALSE]
  }

  out
}
