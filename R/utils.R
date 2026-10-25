# Argument checkers ------------------------------------------------------------

check_whole <- function(x, arg, min = -Inf, call = rlang::caller_env()) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) &&
    x == round(x) && x >= min
  if (!ok) {
    cli::cli_abort(
      "{.arg {arg}} must be a single whole number >= {min}, not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

check_open_unit <- function(x, arg, call = rlang::caller_env()) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0 && x < 1
  if (!ok) {
    cli::cli_abort(
      "{.arg {arg}} must be a single number in (0, 1), not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

check_flag <- function(x, arg, call = rlang::caller_env()) {
  if (!rlang::is_bool(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be {.code TRUE} or {.code FALSE}, not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

describe <- function(x) {
  if (is.numeric(x) && length(x) == 1) {
    format(x)
  } else if (rlang::is_string(x)) {
    encodeString(x, quote = '"')
  } else {
    cli::format_inline("{.obj_type_friendly {x}}")
  }
}

# Data -------------------------------------------------------------------------

# A numeric matrix of finite values from a matrix or a data frame.
as_numeric_matrix <- function(x, arg = "x", call = rlang::caller_env()) {
  if (is.data.frame(x)) {
    is_num <- vapply(x, is.numeric, logical(1))
    if (!all(is_num)) {
      cli::cli_abort(
        "All columns of {.arg {arg}} must be numeric; {.val {names(x)[!is_num]}} {?is/are} not.",
        call = call
      )
    }
    x <- as.matrix(x)
    storage.mode(x) <- "double"
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be a numeric matrix or data frame, not {describe(x)}.",
      call = call
    )
  }
  if (anyNA(x) || any(!is.finite(x))) {
    cli::cli_abort(
      "{.arg {arg}} must not contain missing or infinite values.",
      call = call
    )
  }
  storage.mode(x) <- "double"
  x
}

# The methods of the neighbor search, as used by `knn_cpp()`. The k-d tree is
# faster than the BLAS brute-force search for data with few columns.
use_kd_tree <- function(nn_method, p) {
  switch(nn_method,
    auto = p <= 10,
    kd_tree = TRUE,
    brute_force = FALSE
  )
}

nearest_neighbors <- function(ref, query, k, exclude_self, nn_method = "auto") {
  knn_cpp(ref, query, k, exclude_self, use_kd_tree(nn_method, ncol(ref)))
}

# Leading principal axes of centered data, with the largest loading of each
# axis made positive so that the signs do not depend on the platform.
pca_basis <- function(xc, ndim) {
  if (ncol(xc) <= nrow(xc)) {
    v <- eigen(crossprod(xc), symmetric = TRUE)$vectors[, seq_len(ndim), drop = FALSE]
  } else {
    v <- svd(xc, nu = 0, nv = ndim)$v
  }
  big <- apply(abs(v), 2, which.max)
  flip <- sign(v[cbind(big, seq_len(ndim))])
  flip[flip == 0] <- 1
  sweep(v, 2, flip, "*")
}
