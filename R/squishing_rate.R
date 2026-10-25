#' Squishing rate parameter
#'
#' The rate at which [manifold_sculpting()] scales down the dimensions to
#' remove, for tuning [step_manifold_sculpting()].
#'
#' @param range A two-element vector holding the defaults for the smallest and
#'   largest possible values.
#' @param trans A `trans` object from the `scales` package, or `NULL`.
#'
#' @return A `quant_param` object from the 'dials' package.
#'
#' @examples
#' squishing_rate()
#' @export
squishing_rate <- function(range = c(0.95, 0.995), trans = NULL) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    trans = trans,
    label = c(squishing_rate = "Squishing Rate"),
    finalize = NULL
  )
}
