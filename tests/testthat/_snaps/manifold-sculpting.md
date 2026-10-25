# invalid arguments are errors

    Code
      manifold_sculpting(x, ndim = 3)
    Condition
      Error in `manifold_sculpting()`:
      ! `ndim` must be smaller than the number of columns of `x` (3), not 3.
    Code
      manifold_sculpting(x, ndim = 0)
    Condition
      Error in `manifold_sculpting()`:
      ! `ndim` must be a single whole number >= 1, not 0.
    Code
      manifold_sculpting(x, neighbors = 500)
    Condition
      Error in `manifold_sculpting()`:
      ! `neighbors` must be smaller than the number of rows of `x` (500), not 500.
    Code
      manifold_sculpting(x, squishing_rate = 1)
    Condition
      Error in `manifold_sculpting()`:
      ! `squishing_rate` must be a single number in (0, 1), not 1.
    Code
      manifold_sculpting(x, nstart = 0)
    Condition
      Error in `manifold_sculpting()`:
      ! `nstart` must be a single whole number >= 1, not 0.
    Code
      manifold_sculpting(x, init = "isomap")
    Condition
      Error in `manifold_sculpting()`:
      ! `init` must be one of "pca", not "isomap".
    Code
      manifold_sculpting(x, init = x[1:10, ])
    Condition
      Error in `manifold_sculpting()`:
      ! `init` must have 500 rows and at least 2 columns.
    Code
      manifold_sculpting(x, nn_method = "ball_tree")
    Condition
      Error in `manifold_sculpting()`:
      ! `nn_method` must be one of "auto", "kd_tree", or "brute_force", not "ball_tree".
    Code
      manifold_sculpting(x, verbose = NA)
    Condition
      Error in `manifold_sculpting()`:
      ! `verbose` must be `TRUE` or `FALSE`, not `NA`.
    Code
      manifold_sculpting(letters)
    Condition
      Error in `manifold_sculpting()`:
      ! `x` must be a numeric matrix or data frame, not a character vector.
    Code
      manifold_sculpting(data.frame(a = 1:3, b = letters[1:3]))
    Condition
      Error in `manifold_sculpting()`:
      ! All columns of `x` must be numeric; "b" is not.
    Code
      manifold_sculpting(rbind(x, NA))
    Condition
      Error in `manifold_sculpting()`:
      ! `x` must not contain missing or infinite values.
    Code
      manifold_sculpting(matrix(1, 10, 3), neighbors = 3)
    Condition
      Error in `manifold_sculpting()`:
      ! Every sample of `x` coincides with its nearest neighbors.

