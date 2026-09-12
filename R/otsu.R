#' Otsu statistics for a discrete histogram
#'
#' Computes the between-class variance for every candidate threshold and the
#' total variance. Kept separate from otsu() because binarize() must evaluate
#' separability at a user-supplied threshold, not only at the optimum.
#'
#' Assumes validated input: integer levels in [0, L-1], no NA. Callers
#' (otsu(), binarize()) guarantee this; nothing is re-checked here.
#'
#' @param g Numeric vector or matrix of integer levels in [0, L - 1].
#' @param L Integer. Number of histogram levels.
#'
#' @return A list with \code{sigma_B} (length-L vector of between-class
#'   variances; \code{NaN} at degenerate thresholds where one class is empty)
#'   and \code{sigma_T} (total variance).
#'
#' @noRd
otsu_stats <- function(g, L) {

  i <- seq_len(L)                     # levels 0..L-1 held at indices 1..L
  n <- tabulate(g + 1L, nbins = L)
  p <- n / sum(n)                     # eq. (1)

  w_low  <- cumsum(p)                 # eq. (2)
  w_high <- 1 - w_low                 # eq. (3)
  mu_k   <- cumsum(i * p)             # eq. (7)
  mu_T   <- mu_k[L]                   # eq. (8)

  mu_low  <- mu_k / w_low             # eq. (4)
  # The high class follows by subtraction: its levels are never traversed.
  mu_high <- (mu_T - mu_k) / w_high   # eq. (5)

  list(sigma_B = w_low * w_high * (mu_high - mu_low)^2,     # eq. (14)
       sigma_T = sum((i - mu_T)^2 * p))                     # eq. (15)
}

#' Otsu's optimal threshold
#'
#' Selects the grey level that maximises the between-class variance of the
#' histogram, following Otsu (1979).
#'
#' @param g Numeric vector or matrix of \strong{integer} levels in
#'   \code{[0, L - 1]}. No quantisation is performed here: continuous input
#'   raises an error. Normally called from \code{\link{binarize}}, which owns
#'   the quantisation step.
#' @param L Integer. Number of histogram levels. Defaults to 256.
#'
#' @return A list with:
#'   \describe{
#'     \item{threshold}{Integer. The \strong{last level belonging to the low
#'       class}, so the classes are \code{{0..threshold}} and
#'       \code{{threshold+1 .. L-1}}. To take the high class as foreground,
#'       binarise with \code{g > threshold}, not \code{>=}.
#'       \code{NA_integer_} for degenerate images.}
#'     \item{separability}{Numeric in \code{[0, 1]}. Otsu's criterion
#'       \code{eta* = sigma_B / sigma_T}, comparable across specimens.
#'       \code{NA_real_} for degenerate images.}
#'     \item{status}{\code{"ok"} or \code{"degenerate"}.}
#'   }
#'
#' @details
#' An image is degenerate when \code{sigma_T == 0}, that is when every pixel
#' holds the same value: a fully saturated frame, a black frame, or a corrupt
#' file that the reader accepts without complaint. Such cases return an
#' explicit \code{status} instead of raising an error, so that a batch run can
#' log the offending file and carry on.
#'
#' Ties are broken by \code{which.max}, which returns the first maximum. This
#' differs from \code{EBImage::otsu}, which returns the midpoint of all
#' maxima; the two diverge on histograms with a plateau.
#'
#' @references Otsu, N. (1979). A threshold selection method from gray-level
#'   histograms. \emph{IEEE Transactions on Systems, Man, and Cybernetics}
#'   9(1): 62-66.
#'
#' @seealso \code{\link{binarize}}, \code{\link{to_gray}}
#'
#' @export
otsu <- function(g, L = 256L) {

  if (length(L) != 1L || !is.numeric(L) || is.na(L) || L != round(L) || L < 2L)
    stop("otsu(): 'L' must be a single integer >= 2.", call. = FALSE)
  L <- as.integer(L)

  if (!is.numeric(g))
    stop("otsu(): 'g' must be numeric.", call. = FALSE)

  if (anyNA(g))
    stop("otsu(): 'g' contains NA values.", call. = FALSE)

  rng <- range(g)
  if (rng[1L] < 0 || rng[2L] > L - 1L)
    stop(sprintf("otsu(): 'g' outside [0, %d]: observed [%g, %g].",
                 L - 1L, rng[1L], rng[2L]), call. = FALSE)

  # Without this check tabulate() truncates silently and returns a plausible
  # but wrong threshold. It is the only non-noisy failure mode here, which is
  # why the extra pass over the data is worth paying.
  if (!is.integer(g) && any(g != round(g)))
    stop("otsu(): 'g' must contain integers; quantise before calling.",
         call. = FALSE)

  st <- otsu_stats(g, L)

  if (st$sigma_T == 0 || !any(is.finite(st$sigma_B)))
    return(list(threshold    = NA_integer_,
                separability = NA_real_,
                status       = "degenerate"))

  k <- which.max(st$sigma_B)

  list(threshold    = k - 1L,
       separability = st$sigma_B[k] / st$sigma_T,
       status       = "ok")
}
