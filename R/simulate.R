#' Simulated manifolds
#'
#' The Swiss roll and S-curve manifolds of Gashler, Ventura and Martinez
#' (2008), with their intrinsic (two-dimensional) coordinates, to test
#' dimensionality reduction methods.
#'
#' @param n The number of points.
#' @param noise The standard deviation of Gaussian noise added to the
#'   coordinates `x`, `y` and `z`.
#'
#' @details
#' For the Swiss roll, \eqn{t_i = 8i/n + 2} for \eqn{0 \le i < n}, and
#' \deqn{x_i = t_i \sin(t_i), \quad y_i \sim U(-6, 6), \quad
#' z_i = t_i \cos(t_i).}
#' Its intrinsic coordinates are the arc length along the roll,
#' \eqn{u_i = (\sinh^{-1}(t_i) + t_i \sqrt{t_i^2 + 1}) / 2}, and
#' \eqn{v_i = y_i}.
#'
#' For the S-curve, \eqn{t_i = (2.2i - 0.1)\pi / n}, and
#' \deqn{x_i = t_i, \quad y_i = \sin(t_i), \quad z_i \sim U(0, 2).}
#' Its intrinsic coordinates are the arc length along the curve,
#' \eqn{u_i = \int_0^{t_i} \sqrt{\cos^2(w) + 1} \, dw}, and \eqn{v_i = z_i}.
#' (The paper gives \eqn{v_i = y_i}, a misprint since \eqn{y_i} is a function
#' of \eqn{u_i}.)
#'
#' @return A tibble with `n` rows and the columns `x`, `y`, `z` (the
#'   coordinates in three dimensions) and `u`, `v` (the intrinsic
#'   coordinates). The points are ordered by `u`.
#'
#' @references
#' Gashler, M., Ventura, D. and Martinez, T. (2008). Iterative non-linear
#' dimensionality reduction with manifold sculpting. In *Advances in Neural
#' Information Processing Systems 20*, 513–520.
#'
#' @examples
#' set.seed(1)
#' sim_swiss_roll(5)
#' sim_s_curve(5)
#' @name sim_manifolds
NULL

#' @rdname sim_manifolds
#' @export
sim_swiss_roll <- function(n = 2000, noise = 0) {
  check_whole(n, "n", min = 1)
  check_noise(noise)
  t <- 8 * (seq_len(n) - 1) / n + 2
  y <- stats::runif(n, -6, 6)
  out <- tibble::tibble(
    x = t * sin(t),
    y = y,
    z = t * cos(t),
    u = (asinh(t) + t * sqrt(t^2 + 1)) / 2,
    v = y
  )
  add_noise(out, noise)
}

#' @rdname sim_manifolds
#' @export
sim_s_curve <- function(n = 2000, noise = 0) {
  check_whole(n, "n", min = 1)
  check_noise(noise)
  t <- (2.2 * (seq_len(n) - 1) - 0.1) * pi / n
  z <- stats::runif(n, 0, 2)
  # The arc length is computed on a fine grid, from 0, and interpolated
  arc <- function(w) sqrt(cos(w)^2 + 1)
  grid <- seq(min(t, 0), max(t, 0), length.out = 4001)
  h <- diff(grid)
  cum <- c(0, cumsum(h * (arc(grid[-1]) + arc(grid[-length(grid)])) / 2))
  cum <- cum - stats::approx(grid, cum, xout = 0)$y
  out <- tibble::tibble(
    x = t,
    y = sin(t),
    z = z,
    u = stats::approx(grid, cum, xout = t)$y,
    v = z
  )
  add_noise(out, noise)
}

check_noise <- function(noise, call = rlang::caller_env()) {
  ok <- is.numeric(noise) && length(noise) == 1 && !is.na(noise) &&
    is.finite(noise) && noise >= 0
  if (!ok) {
    cli::cli_abort(
      "{.arg noise} must be a single non-negative number, not {describe(noise)}.",
      call = call
    )
  }
  invisible(noise)
}

add_noise <- function(data, noise) {
  if (noise > 0) {
    n <- nrow(data)
    for (col in c("x", "y", "z")) {
      data[[col]] <- data[[col]] + stats::rnorm(n, sd = noise)
    }
  }
  data
}
