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

* New `vignette("spectroscopy")` explores the NIR spectra of meat of the
  modeldata package with Manifold Sculpting, the tidyverse and tidymodels:
  a two-dimensional embedding shows that the arch of the PCA scores is a
  genuine curve of the spectra, a one-dimensional embedding unfolds it
  into a coordinate along the gradient of fat, and a short comparison with
  principal component regression shows that unfolding does not help to
  predict protein.
