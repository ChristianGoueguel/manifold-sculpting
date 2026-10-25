# columns are matched by name

    Code
      predict(fit, test[c("x", "y")])
    Condition
      Error in `predict()`:
      ! `new_data` is missing the column "z".

# columns are matched by position without names

    Code
      predict(unnamed, new[, 1:2])
    Condition
      Error in `predict()`:
      ! `new_data` must have 3 columns, not 2.

# edge cases

    Code
      predict(fit, test, refine = NA)
    Condition
      Error in `predict()`:
      ! `refine` must be `TRUE` or `FALSE`, not `NA`.
    Code
      predict(fit, test, foo = 1)
    Condition
      Error in `predict()`:
      ! `...` must be empty.
      x Problematic argument:
      * foo = 1

