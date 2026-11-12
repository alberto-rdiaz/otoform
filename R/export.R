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
