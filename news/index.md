# Changelog

## manifoldsculpting 1.2.1

- Initial version.

- [`manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/manifold_sculpting.md)
  implements Manifold Sculpting (Gashler, Ventura and Martinez, 2008),
  following the reference implementation of the authors in the Waffles
  library. The neighbor search and the iterations are written in C++.
  All the dimensions of the data are kept during the iterations at the
  cost of the preserved ones only, and `nstart` runs the algorithm
  several times and keeps the embedding with the lowest error.

- [`predict()`](https://rdrr.io/r/stats/predict.html) embeds new samples
  from their nearest training samples.

- [`step_manifold_sculpting()`](https://christiangoueguel.com/manifold-sculpting/reference/step_manifold_sculpting.md)
  is a recipes step, with the tuning parameters `num_comp`, `neighbors`
  and
  [`squishing_rate()`](https://christiangoueguel.com/manifold-sculpting/reference/squishing_rate.md).

- [`sim_swiss_roll()`](https://christiangoueguel.com/manifold-sculpting/reference/sim_manifolds.md)
  and
  [`sim_s_curve()`](https://christiangoueguel.com/manifold-sculpting/reference/sim_manifolds.md)
  simulate the test manifolds of the paper.

- New
  [`vignette("comparison")`](https://christiangoueguel.com/manifold-sculpting/articles/comparison.md)
  compares Manifold Sculpting with Isomap, LLE, t-SNE, UMAP, nonlinear
  PCA and kernel PCA. It is precomputed, from
  `vignettes/comparison.Rmd.orig`, by `vignettes/precompute.R`.

- New
  [`vignette("spectroscopy")`](https://christiangoueguel.com/manifold-sculpting/articles/spectroscopy.md)
  applies Manifold Sculpting to the NIR spectra of meat of the modeldata
  package, with the tidyverse and tidymodels: it unfolds the arch of the
  PCA scores, and compares principal component regression with
  regression on the embedding for the prediction of protein.
