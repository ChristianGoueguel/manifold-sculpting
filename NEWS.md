# manifoldsculpting 0.1.0

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
