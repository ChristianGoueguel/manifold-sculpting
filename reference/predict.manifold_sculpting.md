# Embed new samples with Manifold Sculpting

Manifold Sculpting has no explicit mapping from the data to the
embedding. A new sample is embedded from its `neighbors` nearest
training samples: it starts at the combination of their embeddings given
by its locally linear reconstruction weights (as in locally linear
embedding), and is then moved to minimize the Manifold Sculpting error,
the discrepancy between its original and embedded distances and angles
to its neighbors, whose embeddings are held fixed. A new sample that is
identical to a training sample takes its embedding, so the training data
are mapped to
[`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md)'s
`embedding`.

## Usage

``` r
# S3 method for class 'manifold_sculpting'
predict(object, new_data, refine = TRUE, ...)
```

## Arguments

- object:

  A
  [`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md)
  fit.

- new_data:

  A numeric matrix or data frame with the variables used to fit
  `object`. Columns are matched by name when `object` was fitted with
  column names, and by position otherwise.

- refine:

  If `FALSE`, the samples stay at their locally linear reconstruction,
  without the error minimization.

- ...:

  Not used.

## Value

A matrix with one row per sample of `new_data` and `object$ndim`
columns.

## Examples

``` r
set.seed(1)
roll <- sim_swiss_roll(1000)
train <- roll[1:800, ]
test <- roll[801:1000, ]
fit <- manifold_sculpting(train[c("x", "y", "z")], nstart = 3)
head(predict(fit, test))
#>            MS1       MS2
#> [1,] -15.65776 13.589965
#> [2,] -20.51788  4.960831
#> [3,] -15.11771 14.832491
#> [4,] -16.49914 12.496426
#> [5,] -19.30644  7.569465
#> [6,] -17.01433 11.848336
```
