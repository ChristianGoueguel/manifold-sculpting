# Squishing rate parameter

The rate at which
[`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md)
scales down the dimensions to remove, for tuning
[`step_manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/step_manifold_sculpting.md).

## Usage

``` r
squishing_rate(range = c(0.95, 0.995), trans = NULL)
```

## Arguments

- range:

  A two-element vector holding the defaults for the smallest and largest
  possible values.

- trans:

  A `trans` object from the `scales` package, or `NULL`.

## Value

A `quant_param` object from the 'dials' package.

## Examples

``` r
squishing_rate()
#> Squishing Rate (quantitative)
#> Range: [0.95, 0.995]
```
