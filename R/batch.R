#' Process a single otolith image
#'
#' Runs the full contour pipeline on one JPEG: luma conversion, binarisation,
#' largest-component selection, boundary tracing, arc-length resampling,
#' smoothing, principal-axis rotation, rostrum anchoring, unit-area
#' normalisation and centroid-distance signature. Any error is caught and
#' reported in the result so that a batch is never interrupted by a single
#' image.
#'
#' @param path Character. Path to a JPEG file.
#' @param method,threshold Binarisation settings, passed unchanged to
#'   \code{\link{binarize}}, which validates them. \code{threshold} is given on
#'   the quantised scale (the last level of the low class), not on
#'   \code{[0, 1]}, and must be \code{NULL} when \code{method = "otsu"}.
#' @param k Integer. Smoothing parameter passed to \code{\link{smoothout}}.
#' @param grid_size Integer. Number of points of the arc-length resampling.
#' @param n_points Integer. Number of points of the final contour.
#' @param min_area_ratio Numeric. Lower bound for the ratio between the area
#'   of the traced polygon and the pixel count of the component it was traced
#'   from. On the 10 reference images the ratio lies in 0.9975-1.007: the
#'   polygon through pixel centres is ~0.25\% smaller than the pixel count, and
#'   interior pores push it slightly above 1. A ratio clearly below 1 can only
#'   mean a partial or wrong polygon, so the image is reported as an error.
#'   The default leaves a margin of about five times the discretisation
#'   effect. Values above 1 are never an error: they measure interior holes in
#'   the mask and are kept in \code{area_ratio} for quality control.
#' @param overlay_file Optional path of a PNG with the contour drawn on the
#'   image, written by \code{\link{write_overlay}}. \code{NULL} writes nothing.
#'
#' @return A list with:
#'   \describe{
#'     \item{file}{Base name of \code{path}.}
#'     \item{raw}{Traced contour in pixel coordinates. \code{NULL} on error.}
#'     \item{points}{Final \code{n_points} x 2 contour in pixels with
#'       attributes \code{status}, \code{s_star}, \code{delta} and
#'       \code{perimeter}. \code{NULL} on error.}
#'     \item{r}{Centroid-distance signature of the contour after
#'       normalisation to unit area: distance from the area centroid to each
#'       of the \code{n_points} vertices, dimensionless. Carries the attribute
#'       \code{centroid}. This is the input of \code{\link{atrous_b3}}.
#'       \code{NULL} on error.}
#'     \item{area}{Area enclosed by \code{points}, in squared pixels, before
#'       normalisation. \code{NA} on error.}
#'     \item{theta}{Rotation angle applied, in radians. \code{NA} on error.}
#'     \item{center}{Area centroid of the smoothed contour, in image pixel
#'       coordinates: the point subtracted by \code{rotation()} and needed to
#'       map \code{points} back onto the image. \code{NA} on error.}
#'     \item{status}{The \code{status} attribute of \code{points}, or
#'       \code{"error"}.}
#'     \item{msg}{Error message; \code{NA} otherwise.}
#'     \item{threshold}{Threshold actually applied, on the quantised scale.
#'       With \code{method = "otsu"} it varies from image to image.}
#'     \item{separability}{Otsu's eta evaluated at the applied threshold.}
#'     \item{inverted}{Whether the mask polarity had to be flipped.}
#'     \item{n_components}{Number of connected components in the mask before
#'       keeping the largest; a quality-control diagnostic.}
#'     \item{area_ratio}{Polygon area over pixel count of the traced
#'       component; see \code{min_area_ratio}.}
#'     \item{method}{Binarisation method used.}
#'     \item{k}{Smoothing parameter used.}
#'   }
#'
#' @seealso \code{\link{process_directory}}, \code{\link{binarize}},
#'   \code{\link{keep_largest}}, \code{\link{Conte}},
#'   \code{\link{write_overlay}}
#' @export
process_image <- function(path, method = c("otsu", "fixed"), threshold = NULL,
                          k = 128L, grid_size = 4096L, n_points = 512L,
                          min_area_ratio = 0.95, overlay_file = NULL) {

  method <- match.arg(method)
  img    <- NULL   # stays NULL if readJPEG fails; then nothing can be drawn

  res <- tryCatch({

    img <- jpeg::readJPEG(path)
    b   <- binarize(to_gray(img), method = method, threshold = threshold)
    if (b$status != "ok")
      stop("binarize(): degenerate image (a single grey level); ",
           "nothing to trace.", call. = FALSE)

    m   <- keep_largest(b$bw)
    raw <- Conte(m)

    area_ratio <- abs(signed_area(raw)) / sum(m)
    if (area_ratio < min_area_ratio)
      stop(sprintf("Conte(): traced polygon covers only %.1f%% of the ",
                   100 * area_ratio),
           "component area; partial or wrong contour.", call. = FALSE)

    fine    <- resample_arc(raw, grid_size)
    smooth  <- smoothout(fine, k)
    rotated <- rotation(smooth)
    points  <- anchor_rostrum(rotated, n_points)

    norm <- normalize_area(points)
    r    <- centroid_distance(norm)

    list(raw          = raw,
         points       = points,
         r            = r,
         area         = attr(norm, "area"),
         theta        = attr(rotated, "theta"),
         center       = attr(rotated, "center"),
         status       = attr(points, "status"),
         msg          = NA_character_,
         threshold    = b$threshold,
         separability = b$separability,
         inverted     = b$inverted,
         n_components = attr(m, "n_components"),
         area_ratio   = area_ratio)

  }, error = function(e) {
    list(raw          = NULL,
         points       = NULL,
         r            = NULL,
         area         = NA_real_,
         theta        = NA_real_,
         center       = c(NA_real_, NA_real_),
         status       = "error",
         msg          = conditionMessage(e),
         threshold    = NA_integer_,
         separability = NA_real_,
         inverted     = NA,
         n_components = NA_integer_,
         area_ratio   = NA_real_)
  })

  # Outside the tryCatch: a failure to write the PNG must not turn a correctly
  # processed otolith into an error.
  if (!is.null(overlay_file) && !is.null(img))
    write_overlay(img, res, overlay_file)

  c(list(file = basename(path)),
    res,
    list(method = method, k = k))
}

#' Process every image in a directory
#'
#' Applies \code{\link{process_image}} to each JPEG in a folder and writes the
#' results next to the images, in \code{out_dir}:
#' \describe{
#'   \item{\code{log.csv}}{One row per image with the scalar fields of the
#'     result (status, threshold, area ratio, ...). Appended after each image,
#'     so an interrupted run keeps what was done.}
#'   \item{\code{results.rds}}{The full list of results, for R.}
#'   \item{\code{overlay/}}{With \code{overlay = TRUE}, one PNG per image with
#'     the contour drawn on it (\code{\link{write_overlay}}).}
#' }
#' Failures are recorded rather than propagated, so one bad image does not
#' stop the run.
#'
#' @param dir Path to the folder containing the images.
#' @param out_dir Where to write the outputs. Defaults to a folder
#'   \code{otoform} inside \code{dir}, so each locality keeps its own results.
#' @param overlay Logical. Write the overlay PNGs?
#' @inheritDotParams process_image -path -overlay_file
#'
#' @return A list of results, one per image, named by file (invisibly: the
#'   same list is in \code{results.rds}).
#' @export
process_directory <- function(dir, out_dir = file.path(dir, "otoform"),
                              overlay = FALSE, ...) {

  if (!dir.exists(dir))
    stop("process_directory(): directory not found: ", dir, call. = FALSE)

  files <- list.files(dir, pattern = "\\.jpe?g$", full.names = TRUE,
                      ignore.case = TRUE)
  if (length(files) == 0L)
    stop("process_directory(): no JPEG files found in ", dir, call. = FALSE)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (overlay)
    dir.create(file.path(out_dir, "overlay"), showWarnings = FALSE)
  log_file <- file.path(out_dir, "log.csv")

  out <- vector("list", length(files))
  names(out) <- basename(files)

  for (i in seq_along(files)) {
    message(sprintf("[%d/%d] %s", i, length(files), basename(files[i])))

    overlay_file <- if (overlay)
      file.path(out_dir, "overlay",
                sub("\\.jpe?g$", ".png", basename(files[i]), ignore.case = TRUE))
    else NULL

    out[[i]] <- process_image(files[i], ..., overlay_file = overlay_file)

    utils::write.table(log_row(out[[i]]), log_file, sep = ",",
                       row.names = FALSE, col.names = i == 1L,
                       append = i > 1L)
  }

  saveRDS(out, file.path(out_dir, "results.rds"))
  invisible(out)
}

#' Wavelet-decompose every signature in a collection
#'
#' Applies \code{\link{atrous_b3}} to the centroid-distance signature
#' \code{r} of each element returned by \code{\link{process_directory}} (or
#' by \code{lapply} over \code{\link{process_image}}). Elements whose
#' contour stage did not succeed are carried through with their status, so
#' the output is always parallel to the input. Failures in the decomposition
#' itself are recorded rather than propagated.
#'
#' @param res List returned by \code{\link{process_directory}}: one element
#'   per image, each carrying \code{file}, \code{r} and \code{status}.
#' @param J Integer. Number of decomposition levels, passed to
#'   \code{\link{atrous_b3}}, which validates it against the signature length.
#'
#' @return A list with one element per input element, named by \code{file}
#'   and carrying attribute \code{J}. Each element is a list with:
#'   \describe{
#'     \item{file}{Image file name.}
#'     \item{c}{Approximations \code{c[[1]]} (the signature itself) to
#'       \code{c[[J + 1]]}. \code{NULL} if not computed.}
#'     \item{w}{Details \code{w[[1]]} to \code{w[[J]]}, with
#'       \code{w[[j]] = c[[j]] - c[[j + 1]]}. \code{NULL} if not computed.}
#'     \item{status}{\code{"ok"}, \code{"error"}, or the status inherited
#'       from the contour stage.}
#'     \item{msg}{Error message; \code{NA} otherwise.}
#'   }
#'
#' @seealso \code{\link{process_directory}}, \code{\link{atrous_b3}}
#' @export
process_contours <- function(res, J = 6L) {

  if (!is.list(res) || length(res) == 0L)
    stop("process_contours(): 'res' must be a non-empty list, as returned ",
         "by process_directory().", call. = FALSE)

  out <- vector("list", length(res))
  names(out) <- vapply(res, function(z) as.character(z$file), character(1L))

  for (i in seq_along(res)) {
    z <- res[[i]]

    if (identical(z$status, "ok") && !is.null(z$r)) {
      out[[i]] <- tryCatch({
        w <- atrous_b3(z$r, J = J)
        list(file = z$file, c = w$c, w = w$w,
             status = "ok", msg = NA_character_)
      }, error = function(e) {
        list(file = z$file, c = NULL, w = NULL,
             status = "error", msg = conditionMessage(e))
      })
    } else {
      out[[i]] <- list(file = z$file, c = NULL, w = NULL,
                       status = if (is.null(z$status)) "error" else z$status,
                       msg    = if (is.null(z$msg)) NA_character_ else z$msg)
    }
  }

  attr(out, "J") <- J
  out
}
