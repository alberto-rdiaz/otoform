#' Process a single otolith image
#'
#' Runs the full pipeline on one image: binarisation, contour extraction,
#' resampling to a uniform grid, smoothing, centring, rotation and
#' sub-pixel rostrum anchoring.
#'
#' Errors are captured rather than propagated, so a failing image does not
#' interrupt batch processing. Check the \code{status} field.
#'
#' @param path Path to a JPEG image.
#' @param threshold Binarisation threshold passed to \code{toBW}.
#' @param k Number of smoothing passes.
#' @param grid_size Size of the intermediate uniform resampling grid.
#' @param n_points Number of points in the final normalised contour.
#'
#' @return A list with the raw contour in image coordinates, the normalised
#'   contour, the rotation angle, the parameters used, and a status field.
#' @export
#'
process_image <- function(path, threshold = 0.5, k = 128L,
                          grid_size = 4096L, n_points = 512L) {

  res <- tryCatch({

    img <- jpeg::readJPEG(path)
    bw  <- toBW(img, threshold)
    raw <- Conte(bw)

    fine    <- resample_arc(raw, grid_size)
    smooth  <- smoothout(fine, k)
    rotated <- rotation(smooth)
    points  <- anchor_rostrum(rotated, n_points)

    list(raw    = raw,
         points = points,
         theta  = attr(rotated, "theta"),
         status = attr(points, "status"),
         msg    = NA_character_)

  }, error = function(e) {
    list(raw = NULL, points = NULL, theta = NA_real_,
         status = "error", msg = conditionMessage(e))
  })

  c(list(file = basename(path)),
    res,
    list(threshold = threshold, k = k))
}

#' Process every image in a directory
#'
#' Applies \code{process_image} to each JPEG in a folder. Failures are
#' recorded rather than propagated, so one bad image does not stop the run.
#'
#' @param dir Path to the folder containing the images.
#' @param ... Arguments passed on to \code{process_image}.
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
