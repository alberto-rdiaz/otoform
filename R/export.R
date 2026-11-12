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
