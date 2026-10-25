#' Manifold Sculpting signal extraction
#'
#' `step_manifold_sculpting()` creates a *specification* of a recipe step that
#' converts numeric data into one or more new dimensions with Manifold
#' Sculpting.
#'
#' @inheritParams recipes::step_pca
#' @param ... One or more selector functions to choose variables for this
#'   step. See [recipes::selections()] for more details.
#' @param num_comp The number of dimensions of the embedding. It is reduced
#'   to one less than the number of selected columns when it is larger. With
#'   `num_comp = 0`, the data are returned unchanged.
#' @param neighbors The number of nearest neighbors that define the local
#'   relationships to preserve.
#' @param squishing_rate The factor by which the dimensions to remove are
#'   scaled down at each iteration, in (0, 1).
#' @param options A list of other arguments to [manifold_sculpting()]:
#'   `nstart`, `max_iter`, `patience`, `nn_method` and `verbose`.
#' @param seed A single integer, the seed of the random numbers used when the
#'   step is trained, for reproducible results.
#' @param res The [manifold_sculpting()] fit, once the step has been trained
#'   by [recipes::prep()].
#' @param prefix The prefix of the names of the new columns, which are
#'   numbered from 1 (`MS1`, `MS2`, ...).
#'
#' @details
#' Manifold Sculpting (Gashler, Ventura and Martinez, 2008) finds a
#' low-dimensional embedding that preserves the distances and angles between
#' neighboring samples; see [manifold_sculpting()] for the algorithm. It is
#' based on Euclidean distances, so the predictors should usually be
#' centered and scaled first, with [recipes::step_normalize()].
#'
#' New samples, given to [recipes::bake()], are embedded from their nearest
#' training samples by [predict.manifold_sculpting()]. The training set gets
#' the embedding found by the algorithm.
#'
#' # Tidying
#'
#' When you [`tidy()`][recipes::tidy.recipe] this step, a tibble is returned
#' with the columns `terms` (the selectors or variables selected) and `id`.
#'
#' # Tuning parameters
#'
#' This step has 3 tuning parameters: `num_comp` (# Components), `neighbors`
#' (# Nearest Neighbors) and `squishing_rate` (Squishing Rate, see
#' [squishing_rate()]).
#'
#' # Case weights
#'
#' The underlying operation does not allow for case weights.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of any existing operations.
#'
#' @references
#' Gashler, M., Ventura, D. and Martinez, T. (2008). Iterative non-linear
#' dimensionality reduction with manifold sculpting. In *Advances in Neural
#' Information Processing Systems 20*, 513–520.
#'
#' @examples
#' library(recipes)
#'
#' set.seed(1)
#' roll <- sim_swiss_roll(1000)
#' train <- roll[1:800, ]
#' test <- roll[801:1000, ]
#'
#' rec <- recipe(~ x + y + z, data = train) |>
#'   step_manifold_sculpting(all_numeric_predictors(), num_comp = 2, seed = 1)
#'
#' rec_trained <- prep(rec)
#' rec_trained
#' bake(rec_trained, new_data = test)
#' tidy(rec_trained, number = 1)
#' @export
step_manifold_sculpting <- function(recipe,
                                    ...,
                                    role = "predictor",
                                    trained = FALSE,
                                    num_comp = 2,
                                    neighbors = 20,
                                    squishing_rate = 0.99,
                                    options = list(),
                                    seed = sample.int(10^5, 1),
                                    res = NULL,
                                    columns = NULL,
                                    prefix = "MS",
                                    keep_original_cols = FALSE,
                                    skip = FALSE,
                                    id = recipes::rand_id("manifold_sculpting")) {
  recipes::add_step(
    recipe,
    step_manifold_sculpting_new(
      terms = rlang::enquos(...),
      role = role,
      trained = trained,
      num_comp = num_comp,
      neighbors = neighbors,
      squishing_rate = squishing_rate,
      options = options,
      seed = seed,
      res = res,
      columns = columns,
      prefix = prefix,
      keep_original_cols = keep_original_cols,
      skip = skip,
      id = id
    )
  )
}

step_manifold_sculpting_new <- function(terms, role, trained, num_comp,
                                        neighbors, squishing_rate, options,
                                        seed, res, columns, prefix,
                                        keep_original_cols, skip, id) {
  recipes::step(
    subclass = "manifold_sculpting",
    terms = terms,
    role = role,
    trained = trained,
    num_comp = num_comp,
    neighbors = neighbors,
    squishing_rate = squishing_rate,
    options = options,
    seed = seed,
    res = res,
    columns = columns,
    prefix = prefix,
    keep_original_cols = keep_original_cols,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_manifold_sculpting <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], types = c("double", "integer"))
  check_whole(x$num_comp, "num_comp", min = 0)
  check_whole(x$neighbors, "neighbors", min = 1)
  check_open_unit(x$squishing_rate, "squishing_rate")
  check_whole(x$seed, "seed")
  check_step_options(x$options)
  if (!rlang::is_string(x$prefix)) {
    cli::cli_abort("{.arg prefix} must be a single string, not {describe(x$prefix)}.")
  }

  num_comp <- x$num_comp
  if (num_comp >= length(col_names) && length(col_names) > 0) {
    num_comp <- length(col_names) - 1
    cli::cli_warn(
      "{.arg num_comp} was reduced to {num_comp}, one less than the number of selected columns."
    )
  }

  if (num_comp > 0 && length(col_names) > 0) {
    args <- c(
      list(
        x = training[, col_names],
        ndim = num_comp,
        neighbors = min(x$neighbors, nrow(training) - 1),
        squishing_rate = x$squishing_rate
      ),
      x$options
    )
    res <- withr::with_seed(x$seed, rlang::exec(manifold_sculpting, !!!args))
  } else {
    res <- NULL
  }

  step_manifold_sculpting_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    num_comp = num_comp,
    neighbors = x$neighbors,
    squishing_rate = x$squishing_rate,
    options = x$options,
    seed = x$seed,
    res = res,
    columns = col_names,
    prefix = x$prefix,
    keep_original_cols = recipes::get_keep_original_cols(x),
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_manifold_sculpting <- function(object, new_data, ...) {
  col_names <- names(object$columns)
  recipes::check_new_data(col_names, object, new_data)
  if (object$num_comp == 0 || length(col_names) == 0) {
    return(new_data)
  }

  comps <- stats::predict(object$res, new_data[, col_names])
  colnames(comps) <- recipes::names0(ncol(comps), object$prefix)
  rownames(comps) <- NULL
  comps <- tibble::as_tibble(comps)
  comps <- recipes::check_name(comps, new_data, object)
  new_data <- vctrs::vec_cbind(new_data, comps, .name_repair = "minimal")
  recipes::remove_original_cols(new_data, object, col_names)
}

#' @export
print.step_manifold_sculpting <- function(x, width = max(20, options()$width - 35), ...) {
  if (x$num_comp == 0) {
    title <- "Manifold Sculpting was not conducted for "
  } else {
    title <- "Manifold Sculpting extraction with "
  }
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @rdname step_manifold_sculpting
#' @param x A `step_manifold_sculpting` object.
#' @export
tidy.step_manifold_sculpting <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble::tibble(terms = unname(x$columns))
  } else {
    res <- tibble::tibble(terms = recipes::sel2char(x$terms))
  }
  res$id <- rep(x$id, nrow(res))
  res
}

#' @export
tunable.step_manifold_sculpting <- function(x, ...) {
  tibble::tibble(
    name = c("num_comp", "neighbors", "squishing_rate"),
    call_info = list(
      list(pkg = "dials", fun = "num_comp", range = c(1L, 10L)),
      list(pkg = "dials", fun = "neighbors", range = c(10L, 40L)),
      list(pkg = "manifoldsculpting", fun = "squishing_rate")
    ),
    source = "recipe",
    component = "step_manifold_sculpting",
    component_id = x$id
  )
}

#' @export
required_pkgs.step_manifold_sculpting <- function(x, ...) {
  "manifoldsculpting"
}

check_step_options <- function(options, call = rlang::caller_env()) {
  allowed <- c("nstart", "max_iter", "patience", "nn_method", "verbose")
  if (!is.list(options) || (length(options) > 0 && !rlang::is_named(options))) {
    cli::cli_abort("{.arg options} must be a named list.", call = call)
  }
  unknown <- setdiff(names(options), allowed)
  if (length(unknown) > 0) {
    cli::cli_abort(
      c(
        "{.arg options} has unknown element{?s} {.val {unknown}}.",
        "i" = "Allowed elements are {.val {allowed}}."
      ),
      call = call
    )
  }
  invisible(options)
}
