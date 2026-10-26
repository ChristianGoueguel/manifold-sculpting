# manifoldsculpting 1.2.1

* Initial version.

* `manifold_sculpting()` implements Manifold Sculpting (Gashler, Ventura and
  Martinez, 2008), following the reference implementation of the authors in
  the Waffles library. The neighbor search and the iterations are written in
  C++. All the dimensions of the data are kept during the iterations at the
  cost of the preserved ones only, and `nstart` runs the algorithm several
  times and keeps the embedding with the lowest error.

* `predict()` embeds new samples from their nearest training samples.

* `step_manifold_sculpting()` is a recipes step, with the tuning parameters
  `num_comp`, `neighbors` and `squishing_rate()`.

* `sim_swiss_roll()` and `sim_s_curve()` simulate the test manifolds of the
  paper.

* New `vignette("comparison")` compares Manifold Sculpting with Isomap, LLE,
  t-SNE, UMAP, nonlinear PCA and kernel PCA. It is precomputed, from
  `vignettes/comparison.Rmd.orig`, by `vignettes/precompute.R`.

* New `vignette("spectroscopy")` applies Manifold Sculpting to the NIR
  spectra of meat of the modeldata package, with the tidyverse and
  tidymodels: it unfolds the arch of the PCA scores, and compares
  principal component regression with regression on the embedding for
  the prediction of protein.
