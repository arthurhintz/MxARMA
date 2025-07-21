#' Random Generation from the Parametrized Maxwell Distribution
#'
#' Generates \code{n} random values from the Maxwell distribution with a given scale parameter \eqn{\mu > 0},
#' using the quantile function.
#'
#' @param n Number of observations to generate. Must be a positive integer.
#' @param mu A positive numeric value representing the scale parameter \eqn{\mu}.
#'
#' @return A numeric vector of random values from the Maxwell distribution.
#'
#' @examples
#' rmax(n = 5, mu = 2)
#'
#' @export
rmax <- function(n, mu) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n))
    stop("'n' must be a positive integer.")
  if (!is.numeric(mu) || length(mu) != 1 || mu <= 0)
    stop("'mu' must be a single positive number.")

  u <- runif(n)
  r <- numeric(n)

  for (i in seq_len(n)) {
    r[i] <- MxARMA::qmax(u[i], mu)
  }
  return(r)
}

