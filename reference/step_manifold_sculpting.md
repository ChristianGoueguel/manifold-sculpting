# Manifold Sculpting signal extraction

`step_manifold_sculpting()` creates a *specification* of a recipe step
that converts numeric data into one or more new dimensions with Manifold
Sculpting.

## Usage

``` r
step_manifold_sculpting(
  recipe,
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
  id = recipes::rand_id("manifold_sculpting")
)

# S3 method for class 'step_manifold_sculpting'
tidy(x, ...)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- ...:

  One or more selector functions to choose variables for this step. See
  [`recipes::selections()`](https://recipes.tidymodels.org/reference/selections.html)
  for more details.

- role:

  For model terms created by this step, what analysis role should they
  be assigned? By default, the new columns created by this step from the
  original variables will be used as *predictors* in a model.

- trained:

  A logical to indicate if the quantities for preprocessing have been
  estimated.

- num_comp:

  The number of dimensions of the embedding. It is reduced to one less
  than the number of selected columns when it is larger. With
  `num_comp = 0`, the data are returned unchanged.

- neighbors:

  The number of nearest neighbors that define the local relationships to
  preserve.

- squishing_rate:

  The factor by which the dimensions to remove are scaled down at each
  iteration, in (0, 1).

- options:

  A list of other arguments to
  [`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md):
  `nstart`, `max_iter`, `patience`, `nn_method` and `verbose`.

- seed:

  A single integer, the seed of the random numbers used when the step is
  trained, for reproducible results.

- res:

  The
  [`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md)
  fit, once the step has been trained by
  [`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html).

- columns:

  A character string of the selected variable names. This field is a
  placeholder and will be populated once
  [`prep()`](https://recipes.tidymodels.org/reference/prep.html) is
  used.

- prefix:

  The prefix of the names of the new columns, which are numbered from 1
  (`MS1`, `MS2`, ...).

- keep_original_cols:

  A logical to keep the original variables in the output. Defaults to
  `FALSE`.

- skip:

  A logical. Should the step be skipped when the recipe is baked by
  [`bake()`](https://recipes.tidymodels.org/reference/bake.html)? While
  all operations are baked when
  [`prep()`](https://recipes.tidymodels.org/reference/prep.html) is run,
  some operations may not be able to be conducted on new data (e.g.
  processing the outcome variable(s)). Care should be taken when using
  `skip = TRUE` as it may affect the computations for subsequent
  operations.

- id:

  A character string that is unique to this step to identify it.

- x:

  A `step_manifold_sculpting` object.

## Value

An updated version of `recipe` with the new step added to the sequence
of any existing operations.

## Details

Manifold Sculpting (Gashler, Ventura and Martinez, 2008) finds a
low-dimensional embedding that preserves the distances and angles
between neighboring samples; see
[`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md)
for the algorithm. It is based on Euclidean distances, so the predictors
should usually be centered and scaled first, with
[`recipes::step_normalize()`](https://recipes.tidymodels.org/reference/step_normalize.html).

New samples, given to
[`recipes::bake()`](https://recipes.tidymodels.org/reference/bake.html),
are embedded from their nearest training samples by
[`predict.manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/predict.manifold_sculpting.md).
The training set gets the embedding found by the algorithm.

## Tidying

When you
[`tidy()`](https://recipes.tidymodels.org/reference/tidy.recipe.html)
this step, a tibble is returned with the columns `terms` (the selectors
or variables selected) and `id`.

## Tuning parameters

This step has 3 tuning parameters: `num_comp` (# Components),
`neighbors` (# Nearest Neighbors) and `squishing_rate` (Squishing Rate,
see
[`squishing_rate()`](https://christiangoueguel.com/manifold-sculpting/reference/squishing_rate.md)).

## Case weights

The underlying operation does not allow for case weights.

## References

Gashler, M., Ventura, D. and Martinez, T. (2008). Iterative non-linear
dimensionality reduction with manifold sculpting. In *Advances in Neural
Information Processing Systems 20*, 513–520.

## Examples

``` r
library(recipes)

set.seed(1)
roll <- sim_swiss_roll(1000)
train <- roll[1:800, ]
test <- roll[801:1000, ]

rec <- recipe(~ x + y + z, data = train) |>
  step_manifold_sculpting(all_numeric_predictors(), num_comp = 2, seed = 1)

rec_trained <- prep(rec)
rec_trained
#> 
#> ── Recipe ──────────────────────────────────────────────────────────────────────
#> 
#> ── Inputs 
#> Number of variables by role
#> predictor: 3
#> 
#> ── Training information 
#> Training data contained 800 data points and no incomplete rows.
#> 
#> ── Operations 
#> • Manifold Sculpting extraction with: x, y, z | Trained
bake(rec_trained, new_data = test)
#> # A tibble: 200 × 2
#>      MS1   MS2
#>    <dbl> <dbl>
#>  1 10.3   18.0
#>  2 17.7   11.5
#>  3  9.35  19.0
#>  4 11.4   17.3
#>  5 15.7   13.5
#>  6 12.1   16.8
#>  7 15.1   14.2
#>  8  9.78  19.1
#>  9 15.3   14.3
#> 10 18.1   11.9
#> # ℹ 190 more rows
tidy(rec_trained, number = 1)
#> # A tibble: 3 × 2
#>   terms id                      
#>   <chr> <chr>                   
#> 1 x     manifold_sculpting_JtE6m
#> 2 y     manifold_sculpting_JtE6m
#> 3 z     manifold_sculpting_JtE6m
```
