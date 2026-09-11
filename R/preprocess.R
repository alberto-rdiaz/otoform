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

#' Convert an RGB image to black and white by luminance thresholding
#'
#' @param pics A 3D array (height x width x 3) with RGB channels between
#'   0 and 1, as returned by \code{jpeg::readJPEG}.
#'
#' @param threshold Numeric between 0 and 1. Luminance values at or above
#'   the threshold become 1 (white); values below become 0 (black).
#'
#' @return A 2D matrix (height x width) of 0s and 1s.
#' @export
toBW <- function(pics, threshold) {
  pics <- 0.2126 * pics[,,1] + 0.7152 * pics[,,2] + 0.0722 * pics[,,3]
  pics[which(pics >= threshold)] <- 1
  pics[which(pics <  threshold)] <- 0
  return(pics)
}

#' Compute the signed area of a closed polygon (shoelace formula)
#'
#' The sign of the result indicates the traversal direction: positive for
#' counter-clockwise, negative for clockwise. Used internally to enforce a
#' consistent contour orientation.
#'
#' @param x,y Numeric vectors of vertex coordinates.
#'
#' @return A single numeric value. Positive if counter-clockwise.
#' @keywords internal
#'
signed_area <- function(x, y) {
  sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y) / 2
}

#' Trace the contour of a binary object (Moore boundary tracing)
#'
#' Follows the outer boundary of a binary shape starting from a point found
#' by scanning inward from the image edge.
#'
#' @param I A 2D binary matrix (0s and 1s), typically the output of
#'   \code{toBW}.
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

  if (signed_area(out[, 1L], out[, 2L]) < 0) {
    out <- out[rev(seq_len(nrow(out))), , drop = FALSE]
  }

  out
}
