# Manifold Sculpting

`manifold_sculpting()` reduces the dimension of numeric data with
Manifold Sculpting (Gashler, Ventura and Martinez, 2008), an iterative
nonlinear dimensionality reduction method that preserves the distances
and angles between neighboring points.

## Usage

``` r
manifold_sculpting(
  x,
  ndim = 2,
  neighbors = 20,
  squishing_rate = 0.99,
  nstart = 1,
  max_iter = 5000,
  patience = 50,
  init = "pca",
  nn_method = c("auto", "kd_tree", "brute_force"),
  verbose = FALSE
)
```

## Arguments

- x:

  A numeric matrix or data frame, with one row per sample. Distances are
  Euclidean, so variables on different scales should be standardized
  first.

- ndim:

  The number of dimensions of the embedding. It must be smaller than the
  number of columns of `x`.

- neighbors:

  The number of nearest neighbors, `k`, that define the local
  relationships to preserve. With too few neighbors, the embedding is
  more likely to fold onto itself; neighborhoods that are too large cut
  across the folds of the manifold.

- squishing_rate:

  The factor, `sigma` in (0, 1), by which the dimensions to remove are
  scaled down at each iteration. Values closer to 1 (0.99 in the paper)
  are slower but follow the manifold more closely.

- nstart:

  The number of runs. The runs start from the same configuration but
  take different random paths, and the embedding with the lowest final
  error is returned. Several runs make the result more reliable on
  manifolds that are hard to unfold (see Details).

- max_iter:

  The maximum number of iterations.

- patience:

  Once the dimensions to remove have shrunk below 1% of their size, the
  iterations stop when the error has not improved for `patience`
  iterations.

- init:

  The starting configuration. `"pca"` (the default) rotates the data
  onto their principal components, as recommended in the paper.
  Alternatively, a numeric matrix with one row per sample and at least
  `ndim` columns, such as the result of another dimensionality reduction
  method: its first `ndim` columns are refined, and its other columns,
  if any, are scaled down.

- nn_method:

  The nearest neighbor search: `"kd_tree"`, `"brute_force"` (computed
  with BLAS), or `"auto"`, which uses a k-d tree for data with up to 10
  columns. Both are exact.

- verbose:

  Print the progress every 50 iterations?

## Value

An object of class `manifold_sculpting`, a list with:

- embedding:

  The embedding, a matrix with one row per sample and `ndim` columns.

- error:

  The total error at each iteration (of the best run).

- iterations:

  The number of iterations (of the best run).

- converged:

  Whether the iterations stopped before `max_iter`.

- final_error:

  The final error of each run.

and the elements needed by
[`predict.manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/predict.manifold_sculpting.md)
to embed new samples.

## Details

The algorithm has five steps:

1.  Find the `k` nearest neighbors of each point.

2.  For each point \\p_i\\ and neighbor \\n\_{ij}\\, record the distance
    \\\delta\_{ij}\\ between them, and the angle \\\theta\_{ij}\\ formed
    by \\p_i\\, \\n\_{ij}\\ and \\m\_{ij}\\, the neighbor of \\n\_{ij}\\
    that is the most colinear with \\p_i\\. The average neighbor
    distance is \\\delta\_{ave}\\.

3.  Rotate the data onto their principal components. The first `ndim`
    components are the preserved dimensions; the others are scaled.

4.  Iterate: multiply the scaled dimensions by `squishing_rate`, then
    scale the preserved dimensions up until the average neighbor
    distance is \\\delta\_{ave}\\ again. Visit the points breadth-first,
    from a random point, and move each one in the preserved dimensions,
    by hill climbing, to minimize the error \$\$\epsilon_i = \sum_j
    w\_{ij} \left\[ \left(\frac{\delta\_{ij} -
    \delta\_{ij}^0}{\delta\_{ave}}\right)^2 + \max\left(0,
    \frac{\theta\_{ij}^0 - \theta\_{ij}}{\pi}\right)^2 \right\]\$\$
    where the superscript 0 marks the original values, and the weights
    \\w\_{ij}\\ are higher for the neighbors that have already been
    moved in this iteration.

5.  Drop the scaled dimensions.

The weighting and the normalization of the error, the adaptive step size
of the hill climbing, and the stopping rule follow the reference
implementation of the authors in the Waffles library. Unlike Waffles,
all the dimensions of the data are kept (as the paper asks): the scaled
dimensions only contribute to the distances and angles through a common
factor, so their contributions are computed once and the iterations only
work on the `ndim` preserved dimensions. Points that the breadth-first
traversal cannot reach are also adjusted at each iteration.

Like any local optimization, the iterations can end in a local minimum,
where part of the manifold is folded over another. The final error tells
these runs apart: it is much larger than for a run that unfolded the
manifold. The risk depends on the data, and decreases with more
`neighbors` and more runs (`nstart`). On 2000 points of the Swiss roll
of
[`sim_swiss_roll()`](https://christiangoueguel.com/manifold-sculpting/reference/sim_manifolds.md),
for example, about one run in seven stays folded with 20 neighbors (more
with fewer points), but the best of three runs is reliable, while the
S-curve of
[`sim_s_curve()`](https://christiangoueguel.com/manifold-sculpting/reference/sim_manifolds.md)
is unfolded with 14 neighbors.

The neighbor search and the iterations are written in C++. The algorithm
uses the random number generator of R; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results.

## References

Gashler, M., Ventura, D. and Martinez, T. (2008). Iterative non-linear
dimensionality reduction with manifold sculpting. In *Advances in Neural
Information Processing Systems 20*, 513–520.
<https://proceedings.neurips.cc/paper/2007/hash/c06d06da9666a219db15cf575aff2824-Abstract.html>

Gashler, M. S. (2011). Waffles: A machine learning toolkit. *Journal of
Machine Learning Research*, 12, 2383–2387.

## See also

[`predict.manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/predict.manifold_sculpting.md)
to embed new samples, and
[`step_manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/step_manifold_sculpting.md)
for the 'recipes' step.

## Examples

``` r
set.seed(1)
curve <- sim_s_curve(800)
fit <- manifold_sculpting(curve[c("x", "y", "z")], ndim = 2)
fit
#> Manifold Sculpting embedding
#> * 800 samples, 3 variables to 2 dimensions
#> * 20 neighbors, squishing rate 0.99
#> * 552 iterations (converged), final error 0.003943

# The embedding recovers the intrinsic coordinates of the S-curve
cor(fit$embedding, curve[c("u", "v")])
#>             u           v
#> MS1 0.9996737 -0.01719187
#> MS2 0.4278502  0.90739189
```
