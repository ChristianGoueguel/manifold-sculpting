#' Manifold Sculpting
#'
#' `manifold_sculpting()` reduces the dimension of numeric data with
#' Manifold Sculpting (Gashler, Ventura and Martinez, 2008), an iterative
#' nonlinear dimensionality reduction method that preserves the distances and
#' angles between neighboring points.
#'
#' @param x A numeric matrix or data frame, with one row per sample. Distances
#'   are Euclidean, so variables on different scales should be standardized
#'   first.
#' @param ndim The number of dimensions of the embedding. It must be smaller
#'   than the number of columns of `x`.
#' @param neighbors The number of nearest neighbors, `k`, that define the
#'   local relationships to preserve. With too few neighbors, the embedding
#'   is more likely to fold onto itself; neighborhoods that are too large cut
#'   across the folds of the manifold.
#' @param nstart The number of runs. The runs start from the same
#'   configuration but take different random paths, and the embedding with
#'   the lowest final error is returned. Several runs make the result more
#'   reliable on manifolds that are hard to unfold (see Details).
#' @param squishing_rate The factor, `sigma` in (0, 1), by which the
#'   dimensions to remove are scaled down at each iteration. Values closer to
#'   1 (0.99 in the paper) are slower but follow the manifold more closely.
#' @param max_iter The maximum number of iterations.
#' @param patience Once the dimensions to remove have shrunk below 1% of
#'   their size, the iterations stop when the error has not improved for
#'   `patience` iterations.
#' @param init The starting configuration. `"pca"` (the default) rotates the
#'   data onto their principal components, as recommended in the paper.
#'   Alternatively, a numeric matrix with one row per sample and at least
#'   `ndim` columns, such as the result of another dimensionality reduction
#'   method: its first `ndim` columns are refined, and its other columns, if
#'   any, are scaled down.
#' @param nn_method The nearest neighbor search: `"kd_tree"`,
#'   `"brute_force"` (computed with BLAS), or `"auto"`, which uses a k-d tree
#'   for data with up to 10 columns. Both are exact.
#' @param verbose Print the progress every 50 iterations?
#'
#' @details
#' The algorithm has five steps:
#'
#' 1. Find the `k` nearest neighbors of each point.
#' 2. For each point \eqn{p_i} and neighbor \eqn{n_{ij}}, record the distance
#'    \eqn{\delta_{ij}} between them, and the angle \eqn{\theta_{ij}} formed
#'    by \eqn{p_i}, \eqn{n_{ij}} and \eqn{m_{ij}}, the neighbor of
#'    \eqn{n_{ij}} that is the most colinear with \eqn{p_i}. The average
#'    neighbor distance is \eqn{\delta_{ave}}.
#' 3. Rotate the data onto their principal components. The first `ndim`
#'    components are the preserved dimensions; the others are scaled.
#' 4. Iterate: multiply the scaled dimensions by `squishing_rate`, then
#'    scale the preserved dimensions up until the average neighbor distance
#'    is \eqn{\delta_{ave}} again. Visit the points breadth-first, from a
#'    random point, and move each one in the preserved dimensions, by hill
#'    climbing, to minimize the error
#'    \deqn{\epsilon_i = \sum_j w_{ij} \left[
#'    \left(\frac{\delta_{ij} - \delta_{ij}^0}{\delta_{ave}}\right)^2 +
#'    \max\left(0, \frac{\theta_{ij}^0 - \theta_{ij}}{\pi}\right)^2
#'    \right]}
#'    where the superscript 0 marks the original values, and the weights
#'    \eqn{w_{ij}} are higher for the neighbors that have already been
#'    moved in this iteration.
#' 5. Drop the scaled dimensions.
#'
#' The weighting and the normalization of the error, the adaptive step size
#' of the hill climbing, and the stopping rule follow the reference
#' implementation of the authors in the Waffles library. Unlike Waffles, all
#' the dimensions of the data are kept (as the paper asks): the scaled
#' dimensions only contribute to the distances and angles through a common
#' factor, so their contributions are computed once and the iterations only
#' work on the `ndim` preserved dimensions. Points that the breadth-first
#' traversal cannot reach are also adjusted at each iteration.
#'
#' Like any local optimization, the iterations can end in a local minimum,
#' where part of the manifold is folded over another. The final error tells
#' these runs apart: it is much larger than for a run that unfolded the
#' manifold. The risk depends on the data, and decreases with more
#' `neighbors` and more runs (`nstart`). On 2000 points of the Swiss roll of
#' [sim_swiss_roll()], for example, about one run in seven stays folded with
#' 20 neighbors (more with fewer points), but the best of three runs is
#' reliable, while the S-curve of [sim_s_curve()] is unfolded with 14
#' neighbors.
#'
#' The neighbor search and the iterations are written in C++. The algorithm
#' uses the random number generator of R; call [set.seed()] first for
#' reproducible results.
#'
#' @return An object of class `manifold_sculpting`, a list with:
#' \describe{
#'   \item{embedding}{The embedding, a matrix with one row per sample and
#'     `ndim` columns.}
#'   \item{error}{The total error at each iteration (of the best run).}
#'   \item{iterations}{The number of iterations (of the best run).}
#'   \item{converged}{Whether the iterations stopped before `max_iter`.}
#'   \item{final_error}{The final error of each run.}
#' }
#' and the elements needed by [predict.manifold_sculpting()] to embed new
#' samples.
#'
#' @references
#' Gashler, M., Ventura, D. and Martinez, T. (2008). Iterative non-linear
#' dimensionality reduction with manifold sculpting. In *Advances in Neural
#' Information Processing Systems 20*, 513–520.
#' <https://proceedings.neurips.cc/paper/2007/hash/c06d06da9666a219db15cf575aff2824-Abstract.html>
#'
#' Gashler, M. S. (2011). Waffles: A machine learning toolkit. *Journal of
#' Machine Learning Research*, 12, 2383–2387.
#'
#' @seealso [predict.manifold_sculpting()] to embed new samples, and
#'   [step_manifold_sculpting()] for the 'recipes' step.
#'
#' @examples
#' set.seed(1)
#' curve <- sim_s_curve(800)
#' fit <- manifold_sculpting(curve[c("x", "y", "z")], ndim = 2)
#' fit
#'
#' # The embedding recovers the intrinsic coordinates of the S-curve
#' cor(fit$embedding, curve[c("u", "v")])
#' @export
manifold_sculpting <- function(x,
                               ndim = 2,
                               neighbors = 20,
                               squishing_rate = 0.99,
                               nstart = 1,
                               max_iter = 5000,
                               patience = 50,
                               init = "pca",
                               nn_method = c("auto", "kd_tree", "brute_force"),
                               verbose = FALSE) {
  x <- as_numeric_matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  check_whole(ndim, "ndim", min = 1)
  if (ndim >= p) {
    cli::cli_abort(
      "{.arg ndim} must be smaller than the number of columns of {.arg x} ({p}), not {ndim}."
    )
  }
  check_whole(neighbors, "neighbors", min = 1)
  if (neighbors >= n) {
    cli::cli_abort(
      "{.arg neighbors} must be smaller than the number of rows of {.arg x} ({n}), not {neighbors}."
    )
  }
  check_open_unit(squishing_rate, "squishing_rate")
  check_whole(nstart, "nstart", min = 1)
  check_whole(max_iter, "max_iter", min = 1)
  check_whole(patience, "patience", min = 1)
  nn_method <- rlang::arg_match(nn_method)
  check_flag(verbose, "verbose")

  center <- colMeans(x)
  xc <- x - rep(center, each = n)
  nn <- nearest_neighbors(xc, xc, neighbors, exclude_self = TRUE, nn_method = nn_method)
  if (all(nn$distance == 0)) {
    cli::cli_abort("Every sample of {.arg x} coincides with its nearest neighbors.")
  }

  if (is.character(init)) {
    init <- rlang::arg_match(init, "pca")
    rotation <- pca_basis(xc, ndim)
    y0 <- xc %*% rotation
    resid <- xc - tcrossprod(y0, rotation)
  } else {
    init <- as_numeric_matrix(init, arg = "init")
    if (nrow(init) != n || ncol(init) < ndim) {
      cli::cli_abort(
        "{.arg init} must have {n} rows and at least {ndim} column{?s}."
      )
    }
    rotation <- NULL
    y0 <- init[, seq_len(ndim), drop = FALSE]
    resid <- init[, -seq_len(ndim), drop = FALSE]
  }

  res <- ms_fit_cpp(
    xc, y0, resid, nn$index, nn$distance,
    squishing_rate = squishing_rate,
    max_iter = as.integer(max_iter),
    patience = as.integer(patience),
    burn_in_scale = 0.01,
    nstart = as.integer(nstart),
    verbose = verbose
  )

  embedding <- res$embedding
  colnames(embedding) <- recipes::names0(ndim, "MS")
  rownames(embedding) <- rownames(x)

  structure(
    list(
      embedding = embedding,
      error = res$error,
      iterations = res$iterations,
      converged = res$converged,
      final_error = res$final_error,
      ndim = as.integer(ndim),
      neighbors = as.integer(neighbors),
      squishing_rate = squishing_rate,
      nn_method = nn_method,
      center = center,
      rotation = rotation,
      x = unname(xc),
      nn_index = nn$index,
      nn_distance = nn$distance,
      ave_dist = res$ave_dist
    ),
    class = "manifold_sculpting"
  )
}

#' @export
print.manifold_sculpting <- function(x, ...) {
  status <- if (x$converged) "converged" else "stopped at max_iter"
  cat("Manifold Sculpting embedding\n")
  cat(
    "* ", nrow(x$embedding), " samples, ", length(x$center), " variables to ",
    x$ndim, " dimension", if (x$ndim > 1) "s", "\n",
    "* ", x$neighbors, " neighbors, squishing rate ", x$squishing_rate, "\n",
    "* ", x$iterations, " iterations (", status, "), final error ",
    signif(x$error[length(x$error)], 4), "\n",
    sep = ""
  )
  if (length(x$final_error) > 1) {
    cat(
      "* Best of ", length(x$final_error), " runs, with final errors ",
      paste(signif(x$final_error, 4), collapse = ", "), "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Embed new samples with Manifold Sculpting
#'
#' Manifold Sculpting has no explicit mapping from the data to the
#' embedding. A new sample is embedded from its `neighbors` nearest training
#' samples: it starts at the combination of their embeddings given by its
#' locally linear reconstruction weights (as in locally linear embedding),
#' and is then moved to minimize the Manifold Sculpting error, the
#' discrepancy between its original and embedded distances and angles to
#' its neighbors, whose embeddings are held fixed. A new sample that is
#' identical to a training sample takes its embedding, so the training data
#' are mapped to [manifold_sculpting()]'s `embedding`.
#'
#' @param object A [manifold_sculpting()] fit.
#' @param new_data A numeric matrix or data frame with the variables used to
#'   fit `object`. Columns are matched by name when `object` was fitted with
#'   column names, and by position otherwise.
#' @param refine If `FALSE`, the samples stay at their locally linear
#'   reconstruction, without the error minimization.
#' @param ... Not used.
#'
#' @return A matrix with one row per sample of `new_data` and `object$ndim`
#'   columns.
#'
#' @examples
#' set.seed(1)
#' roll <- sim_swiss_roll(1000)
#' train <- roll[1:800, ]
#' test <- roll[801:1000, ]
#' fit <- manifold_sculpting(train[c("x", "y", "z")], nstart = 3)
#' head(predict(fit, test))
#' @export
predict.manifold_sculpting <- function(object, new_data, refine = TRUE, ...) {
  rlang::check_dots_empty()
  check_flag(refine, "refine")
  x <- new_data_matrix(object, new_data)
  if (nrow(x) == 0) {
    return(object$embedding[0, , drop = FALSE])
  }
  xc <- x - rep(object$center, each = nrow(x))
  nn <- nearest_neighbors(
    object$x, xc, object$neighbors,
    exclude_self = FALSE, nn_method = object$nn_method
  )
  out <- ms_project_cpp(
    object$x, object$embedding, object$nn_index, object$nn_distance,
    xc, nn$index, nn$distance,
    ave_dist = object$ave_dist, refine = refine
  )
  colnames(out) <- colnames(object$embedding)
  rownames(out) <- rownames(x)
  out
}

new_data_matrix <- function(object, new_data, call = rlang::caller_env()) {
  vars <- names(object$center)
  if (!is.null(vars) && !is.null(colnames(new_data))) {
    missing <- setdiff(vars, colnames(new_data))
    if (length(missing) > 0) {
      cli::cli_abort(
        "{.arg new_data} is missing the column{?s} {.val {missing}}.",
        call = call
      )
    }
    new_data <- new_data[, vars, drop = FALSE]
  }
  x <- as_numeric_matrix(new_data, arg = "new_data", call = call)
  if (ncol(x) != length(object$center)) {
    cli::cli_abort(
      "{.arg new_data} must have {length(object$center)} column{?s}, not {ncol(x)}.",
      call = call
    )
  }
  x
}
