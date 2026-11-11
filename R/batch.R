#' Process a single otolith image
#'
#' Runs the full contour pipeline on one JPEG: luma conversion, binarisation,
#' boundary tracing, arc-length resampling, smoothing, principal-axis rotation
#' and rostrum anchoring. Any error is caught and reported in the result so
#' that a batch is never interrupted by a single image.
#'
#' @param path Character. Path to a JPEG file.
#' @param method,threshold Binarisation settings, passed unchanged to
#'   \code{\link{binarize}}, which validates them. \code{threshold} is given on
#'   the quantised scale (the last level of the low class), not on
#'   \code{[0, 1]}, and must be \code{NULL} when \code{method = "otsu"}.
#' @param k Integer. Smoothing parameter passed to \code{\link{smoothout}}.
#' @param grid_size Integer. Number of points of the arc-length resampling.
#' @param n_points Integer. Number of points of the final contour.
#'
#' @return A list with:
#'   \describe{
#'     \item{file}{Base name of \code{path}.}
#'     \item{raw}{Traced contour in pixel coordinates. \code{NULL} on error.}
#'     \item{points}{Final \code{n_points} x 2 contour with attributes
#'       \code{status}, \code{s_star}, \code{delta} and \code{perimeter}.
#'       \code{NULL} on error.}
#'     \item{k}{Smoothing parameter used.}
#'     \item{r}{Centroid-distance signature of the contour after
#'       normalisation to unit area: distance from the area centroid to each
#'       of the \code{n_points} vertices, dimensionless. Carries the attribute
#'       \code{centroid}. This is the input of \code{\link{atrous_b3}}.
#'       \code{NULL} on error.}
#'     \item{area}{Area enclosed by \code{points}, in squared pixels, before
#'       normalisation. \code{NA} on error.}
#'     \item{theta}{Rotation angle applied, in radians. \code{NA} on error.}
#'     \item{status}{The \code{status} attribute of \code{points}, or
#'       \code{"error"}.}
#'     \item{msg}{Error message; \code{NA} otherwise.}
#'     \item{threshold}{Threshold actually applied, on the quantised scale.
#'       With \code{method = "otsu"} it varies from image to image.}
#'     \item{separability}{Otsu's eta evaluated at the applied threshold.}
#'     \item{inverted}{Whether the mask polarity had to be flipped.}
#'     \item{method}{Binarisation method used.}
#'   }
#'
#' @seealso \code{\link{process_directory}}, \code{\link{binarize}}
#' @export
process_image <- function(path, method = c("otsu", "fixed"), threshold = NULL,
                          k = 128L, grid_size = 4096L, n_points = 512L) {

  method <- match.arg(method)

  res <- tryCatch({

    img <- jpeg::readJPEG(path)
    b   <- binarize(to_gray(img), method = method, threshold = threshold)
    if (b$status != "ok")
      stop("binarize(): degenerate image (a single grey level); ",
           "nothing to trace.", call. = FALSE)
    raw <- Conte(b$bw)

    fine    <- resample_arc(raw, grid_size)
    smooth  <- smoothout(fine, k)
    rotated <- rotation(smooth)
    points  <- anchor_rostrum(rotated, n_points)
    norm    <- normalize_area(points)
    r       <- centroid_distance(norm)

    list(raw          = raw,
         points       = points,
         theta        = attr(rotated, "theta"),
         status       = attr(points, "status"),
         msg          = NA_character_,
         threshold    = b$threshold,
         separability = b$separability,
         inverted     = b$inverted,
         r            = r,
         area         = attr(norm, "area"))

  }, error = function(e) {
    list(raw          = NULL,
         points       = NULL,
         theta        = NA_real_,
         status       = "error",
         msg          = conditionMessage(e),
         threshold    = NA_integer_,
         separability = NA_real_,
         inverted     = NA,
         r            = NULL,
         area         = NA_real_)
  })

  c(list(file = basename(path)),
    res,
    list(method = method, k = k))
}

#' Process every image in a directory
#'
#' Applies \code{process_image} to each JPEG in a folder. Failures are
#' recorded rather than propagated, so one bad image does not stop the run.
#'
#' @param dir Path to the folder containing the images.
#' @inheritDotParams process_image -path
#'
#' @return A list of results, one per image, named by file.
#' @export
#'
process_directory <- function(dir, ...) {
  stopifnot(dir.exists(dir))

  files <- list.files(dir, pattern = "\\.jpe?g$",
                      full.names = TRUE, ignore.case = TRUE)

  if (length(files) == 0L) stop("no JPEG files found in ", dir)

  out <- vector("list", length(files))
  names(out) <- basename(files)

  for (i in seq_along(files)) {
    message(sprintf("[%d/%d] %s", i, length(files), basename(files[i])))
    out[[i]] <- process_image(files[i], ...)
  }

  out
}


