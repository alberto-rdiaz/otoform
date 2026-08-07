#' Convert an RGB image to black and white by luminance thresholding
#'
#' @param pics A 3D array (height x width x 3) with RGB channels in [0, 1],
#'   as returned by \code{jpeg::readJPEG}.
#' @param threshold Numeric in [0, 1]. Luminance values at or above the
#'   threshold become 1 (white); values below become 0 (black).
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
#' @return A list with elements X and Y: the pixel coordinates of the
#'   traced contour, in counter-clockwise order.
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

  out <- list(X = Y[-1], Y = (dim(I)[1] - X)[-1])

  if (signed_area(out$X, out$Y) < 0) {
    out$X <- rev(out$X)
    out$Y <- rev(out$Y)
  }

  out
}
