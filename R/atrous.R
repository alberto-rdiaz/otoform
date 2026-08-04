#' À trous wavelet decomposition with B3-spline kernel
#'
#' @param x Numeric vector. The signal to decompose.
#' @param J Integer. Number of decomposition levels. Safe maximum is J = 7
#'   for N = 512. Higher values trigger an error unless allow_unsafe = TRUE.
#' @param allow_unsafe Logical. If FALSE (default), aborts when J is large
#'   enough to cause aliasing. Set TRUE to compute unsafe levels deliberately.
#'
#' @return A list with elements c (approximations) and w (detail planes).
#' @export
atrous_b3 <- function(x, J = 7, allow_unsafe = FALSE) {
  x <- as.numeric(x)
  N <- length(x)

  if (2^J >= N/2 && !allow_unsafe) {
    stop("J=", J, " produce aliasing con N=", N,
         ". Usa allow_unsafe=TRUE si es intencional.")
  }

  h <- c(1, 4, 6, 4, 1) / 16

  periodic_convolve <- function(sig, h, step) {
    L <- length(sig); y <- numeric(L); half <- (length(h)-1)/2
    for (n in 1:L) {
      acc <- 0
      for (k in 1:length(h)) {
        idx <- ((n - 1) + (k - (half+1))*step) %% L + 1
        acc <- acc + h[k] * sig[idx]
      }
      y[n] <- acc
    }
    y
  }

  c_list <- list(x)
  for (j in 1:J) {
    step <- 2^(j-1)
    c_list[[j+1]] <- periodic_convolve(c_list[[j]], h, step)
  }
  w_list <- lapply(1:J, function(j) c_list[[j]] - c_list[[j+1]])

  list(c = c_list, w = w_list)
}
