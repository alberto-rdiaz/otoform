# ---- helpers (internal) ----------------------------------------------------

ok_only <- function(res) Filter(function(z) identical(z$status, "ok"), res)

#' Wide table: one row per element, one column per value
#' @param ok List of successful results.
#' @param get Function extracting a numeric vector of fixed length from one.
#' @param prefix Column-name prefix; columns are \code{prefix_001}, ...
#' @keywords internal
wide_table <- function(ok, get, prefix) {
  m <- do.call(rbind, lapply(ok, function(z) as.numeric(get(z))))
  colnames(m) <- sprintf("%s_%03d", prefix, seq_len(ncol(m)))
  data.frame(file = vapply(ok, `[[`, character(1L), "file"), m,
             check.names = FALSE, stringsAsFactors = FALSE)
}

#' One log row from a process_image() result
#' @param res One element of the list returned by
#'   \code{\link{process_directory}}.
#' @return A one-row data frame with the scalar fields.
#' @keywords internal
log_row <- function(res) {
  data.frame(file = res$file, status = res$status, msg = res$msg,
             method = res$method, threshold = res$threshold,
             separability = res$separability, inverted = res$inverted,
             n_components = res$n_components, area_ratio = res$area_ratio,
             area = res$area, theta = res$theta, k = res$k,
             stringsAsFactors = FALSE)
}

# ---- export_* ---------------------------------------------------------------

#' Export the processing log
#'
#' One row per image, successful or not, with the scalar fields of
#' \code{\link{process_image}}: status, message, binarisation settings and
#' quality-control diagnostics.
#'
#' @param res List returned by \code{\link{process_directory}}.
#' @param file Path of the CSV to write.
#' @return The table, invisibly.
#' @export
export_log <- function(res, file) {
  tab <- do.call(rbind, lapply(res, log_row))
  utils::write.csv(tab, file, row.names = FALSE)
  invisible(tab)
}

#' Export the aligned contours
#'
#' Long table, one row per vertex: \code{file}, \code{i} (vertex index,
#' rostrum at 1), \code{x}, \code{y}. Coordinates are those of
#' \code{points}: pixel units, centred on the area centroid, major axis on x,
#' rostrum on the positive side. This is the aligned contour used for shape
#' comparison, not the position in the source image; \code{theta} and
#' \code{center} in the log map it back (\code{\link{to_image}}). Images with
#' \code{status != "ok"} are left out.
#'
#' @param res List returned by \code{\link{process_directory}}.
#' @param file Path of the CSV to write.
#' @param digits Decimal places kept; sub-pixel anchoring is good to about
#'   0.01 px, so 4 is already more than the data carry.
#' @return The table, invisibly; \code{NULL} if no image succeeded.
#' @export
export_contours <- function(res, file, digits = 4L) {
  ok <- ok_only(res)
  if (length(ok) == 0L) {
    message("export_contours(): no successful image; nothing written.")
    return(invisible(NULL))
  }
  tab <- do.call(rbind, lapply(ok, function(z)
    data.frame(file = z$file, i = seq_len(nrow(z$points)),
               x = round(z$points[, 1L], digits),
               y = round(z$points[, 2L], digits),
               stringsAsFactors = FALSE)))
  utils::write.csv(tab, file, row.names = FALSE)
  invisible(tab)
}

#' Export the centroid-distance signatures
#'
#' Wide table, one row per otolith: \code{file}, \code{r_001} ... \code{r_N}.
#' Dimensionless (unit area), referenced to the area centroid, rostrum at
#' \code{r_001}. The input of \code{\link{atrous_b3}}; full precision.
#'
#' @inheritParams export_contours
#' @return The table, invisibly; \code{NULL} if no image succeeded.
#' @export
export_signatures <- function(res, file) {
  ok <- ok_only(res)
  if (length(ok) == 0L) {
    message("export_signatures(): no successful image; nothing written.")
    return(invisible(NULL))
  }
  tab <- wide_table(ok, function(z) z$r, "r")
  utils::write.csv(tab, file, row.names = FALSE)
  invisible(tab)
}

#' Export the wavelet decomposition, one file per level
#'
#' Writes \code{<prefix>_cJ.csv} (coarse approximation) and
#' \code{<prefix>_w1.csv} ... \code{<prefix>_wJ.csv} (details) into
#' \code{dir}. Each is a wide table, one row per otolith, columns
#' \code{cJ_001...} or \code{wj_001...}; full precision. The signature is
#' recovered exactly as \code{cJ + w1 + ... + wJ}. Elements with
#' \code{status != "ok"} are left out.
#'
#' @param dec List returned by \code{\link{process_contours}}.
#' @param dir Folder where the files are written.
#' @param prefix File-name prefix.
#' @return The paths written, invisibly; \code{NULL} if nothing succeeded.
#' @export
export_wavelets <- function(dec, dir, prefix = "wavelet") {
  J  <- attr(dec, "J")
  ok <- ok_only(dec)
  if (length(ok) == 0L) {
    message("export_wavelets(): no successful element; nothing written.")
    return(invisible(NULL))
  }
  paths <- character(0L)

  f   <- file.path(dir, sprintf("%s_c%d.csv", prefix, J))
  tab <- wide_table(ok, function(z) z$c[[J + 1L]], sprintf("c%d", J))
  utils::write.csv(tab, f, row.names = FALSE)
  paths <- c(paths, f)

  for (j in seq_len(J)) {
    f   <- file.path(dir, sprintf("%s_w%d.csv", prefix, j))
    tab <- wide_table(ok, function(z) z$w[[j]], sprintf("w%d", j))
    utils::write.csv(tab, f, row.names = FALSE)
    paths <- c(paths, f)
  }
  invisible(paths)
}

#' Export the image with the extracted contour drawn on it
#'
#' Writes a PNG of the image with the final contour \code{points} mapped back
#' to image coordinates (red dots), the raw traced contour (grey line), the
#' rostrum anchor (first vertex, green) and the area centroid (black). With
#' \code{status != "ok"} only the image and, if available, \code{raw} are
#' drawn.
#'
#' @param img Image array as returned by \code{jpeg::readJPEG}.
#' @param res One element of the list returned by
#'   \code{\link{process_directory}}.
#' @param file Path of the PNG to write.
#' @param width Width of the PNG in pixels; height follows the image aspect.
#' @return \code{file}, invisibly.
#' @export
export_overlay <- function(img, res, file, width = 1600L) {
  h <- dim(img)[1L]; w <- dim(img)[2L]
  grDevices::png(file, width = width, height = round(width * h / w))
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(0, 0, 0, 0))
  plot(grDevices::as.raster(img))

  if (!is.null(res$raw))
    graphics::lines(res$raw[, 1L] - 0.5, res$raw[, 2L] + 0.5, col = "grey70", lwd = 1)
  if (identical(res$status, "ok")) {
    p <- to_image(res$points, res$theta, res$center)
    graphics::points(p[, 1L] - 0.5, p[, 2L] + 0.5, pch = 19, col = "purple", cex = 1)
    graphics::points(p[1L, 1L] - 0.5, p[1L, 2L] + 0.5, pch = 19, col = "blue", cex = 1) # rostrum
    graphics::points(res$center[1L] - 0.5, res$center[2L] + 0.5, pch = 19, col = "black", cex = 1) # centroid                   # centroid
  }
  invisible(file)
}
