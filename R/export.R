#' Write processing results to disk
#'
#' Writes three files: the complete results as RDS, the normalised
#' coordinates as a long-format text table, and a processing log with one
#' row per image.
#'
#' @param res Output of \code{process_directory}.
#' @param dir Output directory. Created if it does not exist.
#' @param prefix Prefix for the output file names.
#'
#' @return Invisibly, the paths of the three files written.
#' @export
#'
write_results <- function(res, dir, prefix = "otoform") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  f_rds <- file.path(dir, paste0(prefix, "_results.rds"))
  f_txt <- file.path(dir, paste0(prefix, "_coords.txt"))
  f_log <- file.path(dir, paste0(prefix, "_log.txt"))

  saveRDS(res, f_rds)
  utils::write.table(coords_long(res), f_txt,
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(build_log(res), f_log,
                     sep = "\t", row.names = FALSE, quote = FALSE)

  invisible(c(rds = f_rds, txt = f_txt, log = f_log))
}


#' Build a long-format table of normalised coordinates
#'
#' @param res Output of \code{process_directory}.
#'
#' @return A data frame with columns file, index, x and y. Images that
#'   failed are omitted.
#' @export
#'
coords_long <- function(res) {
  ok <- Filter(function(x) !is.null(x$points), res)
  if (length(ok) == 0L) stop("no successful results to export")

  do.call(rbind, lapply(ok, function(x) {
    data.frame(file  = x$file,
               index = seq_len(nrow(x$points)),
               x     = x$points[, 1L],
               y     = x$points[, 2L],
               stringsAsFactors = FALSE)
  }))
}

#' Build a processing log
#'
#' @param res Output of \code{process_directory}.
#'
#' @return A data frame with one row per image: status, parameters used,
#'   number of raw contour vertices, and rotation angle.
#' @export
#'
build_log <- function(res) {
  do.call(rbind, lapply(res, function(x) {
    data.frame(file      = x$file,
               status    = x$status,
               msg       = x$msg,
               n_raw     = if (is.null(x$raw)) NA_integer_ else nrow(x$raw),
               theta_deg = round(x$theta * 180 / pi, 3),
               threshold = x$threshold,
               k         = x$k,
               stringsAsFactors = FALSE)
  }))
}

#' One log row from a process_image() result
#'
#' Keeps the scalar fields of a result (everything but \code{raw},
#' \code{points} and \code{r}) as a one-row data frame, ready to be appended
#' to a CSV or bound into a log.
#'
#' @param res One element of the list returned by
#'   \code{\link{process_directory}}.
#' @return A one-row data frame.
#' @keywords internal
log_row <- function(res) {
  data.frame(file         = res$file,
             status       = res$status,
             msg          = res$msg,
             method       = res$method,
             threshold    = res$threshold,
             separability = res$separability,
             inverted     = res$inverted,
             n_components = res$n_components,
             area_ratio   = res$area_ratio,
             area         = res$area,
             theta        = res$theta,
             k            = res$k,
             stringsAsFactors = FALSE)
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

#' Draw the extracted contour on the source image
#'
#' Writes a PNG of the image with the final contour \code{points} mapped back
#' to image coordinates (red), the raw traced contour (grey), the rostrum
#' anchor (first vertex of \code{points}, green) and the area centroid (blue).
#' With \code{status != "ok"} only the image and, if available, \code{raw}
#' are drawn.
#'
#' @param img Image array as returned by \code{jpeg::readJPEG}.
#' @param res One element of the list returned by
#'   \code{\link{process_directory}}.
#' @param out_file Path of the PNG to write.
#' @param width Width of the PNG in pixels; height follows the image aspect.
#' @return \code{out_file}, invisibly.
#' @export
write_overlay <- function(img, res, out_file, width = 1600L) {
  h <- dim(img)[1L]
  w <- dim(img)[2L]

  grDevices::png(out_file, width = width, height = round(width * h / w))
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(0, 0, 0, 0))
  plot(grDevices::as.raster(img))

  if (!is.null(res$raw))
    graphics::lines(res$raw[, 1L] - 0.5, res$raw[, 2L] + 0.5,
                    col = "grey70", lwd = 1)
  if (identical(res$status, "ok")) {
    p <- to_image(res$points, res$theta, res$center)
    graphics::points(p[, 1L] - 0.5, p[, 2L] + 0.5, pch = 19, col = "purple", cex = 1)
    graphics::points(p[1L, 1L] - 0.5, p[1L, 2L] + 0.5, pch = 19, col = "blue", cex = 1) # rostrum
    graphics::points(res$center[1L] - 0.5, res$center[2L] + 0.5, pch = 19, col = "black", cex = 1) # centroid
  }
  invisible(out_file)
}
