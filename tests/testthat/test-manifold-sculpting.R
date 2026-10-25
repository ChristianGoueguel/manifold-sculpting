curve <- local({
  set.seed(10)
  sim_s_curve(500)
})
x <- as.matrix(curve[c("x", "y", "z")])
fit <- local({
  set.seed(1)
  manifold_sculpting(x, ndim = 2, neighbors = 14)
})

test_that("the embedding recovers the intrinsic coordinates", {
  expect_s3_class(fit, "manifold_sculpting")
  expect_identical(dim(fit$embedding), c(500L, 2L))
  expect_identical(colnames(fit$embedding), c("MS1", "MS2"))
  expect_true(fit$converged)
  expect_gt(procrustes_r2(fit$embedding, curve[c("u", "v")]), 0.995)
  # Manifold Sculpting preserves the distances, so no rescaling is needed
  expect_equal(sum(apply(fit$embedding, 2, var)), var(curve$u) + var(curve$v),
               tolerance = 0.02)
  expect_equal(colMeans(fit$embedding), c(MS1 = 0, MS2 = 0), tolerance = 1e-8)
  expect_length(fit$error, fit$iterations)
})

test_that("it does better than PCA on the Swiss roll", {
  set.seed(3)
  roll <- sim_swiss_roll(800)
  xr <- as.matrix(roll[c("x", "y", "z")])
  set.seed(4)
  ms <- manifold_sculpting(xr, nstart = 2)
  pca <- prcomp(xr)$x[, 1:2]
  expect_gt(procrustes_r2(ms$embedding, roll[c("u", "v")]), 0.98)
  expect_lt(procrustes_r2(pca, roll[c("u", "v")]), 0.5)
})

test_that("results are reproducible with set.seed()", {
  set.seed(1)
  again <- manifold_sculpting(x, ndim = 2, neighbors = 14)
  expect_identical(again$embedding, fit$embedding)
})

test_that("data frames, matrices and both neighbor searches give the same fit", {
  set.seed(1)
  from_df <- manifold_sculpting(curve[c("x", "y", "z")], ndim = 2, neighbors = 14)
  expect_identical(from_df$embedding, fit$embedding)
  set.seed(1)
  brute <- manifold_sculpting(x, ndim = 2, neighbors = 14, nn_method = "brute_force")
  expect_equal(brute$embedding, fit$embedding)
})

test_that("`nstart` keeps the run with the lowest final error", {
  set.seed(5)
  multi <- manifold_sculpting(x, neighbors = 14, nstart = 3, max_iter = 200)
  expect_length(multi$final_error, 3)
  expect_identical(multi$error[multi$iterations], min(multi$final_error))
  expect_output(print(multi), "Best of 3 runs")
})

test_that("`max_iter` stops the iterations", {
  set.seed(1)
  short <- manifold_sculpting(x, neighbors = 14, max_iter = 20)
  expect_identical(short$iterations, 20L)
  expect_false(short$converged)
  expect_output(print(short), "stopped at max_iter")
})

test_that("a curve is unrolled to one dimension", {
  t <- seq(0, 1.5 * pi, length.out = 300)
  spiral <- cbind(t * cos(t), t * sin(t))
  set.seed(1)
  one <- manifold_sculpting(spiral, ndim = 1, neighbors = 10)
  expect_identical(dim(one$embedding), c(300L, 1L))
  expect_gt(abs(cor(one$embedding[, 1], t, method = "spearman")), 0.999)
})

test_that("another embedding can be refined with `init`", {
  set.seed(1)
  start <- as.matrix(curve[c("u", "v")]) + matrix(rnorm(1000, sd = 0.05), ncol = 2)
  refined <- manifold_sculpting(x, neighbors = 14, init = start, max_iter = 100)
  expect_gt(procrustes_r2(refined$embedding, curve[c("u", "v")]), 0.995)
  expect_null(refined$rotation)

  # Extra columns are the dimensions to scale down
  set.seed(1)
  pcs <- prcomp(x)$x
  from_pcs <- manifold_sculpting(x, neighbors = 14, init = pcs)
  expect_gt(procrustes_r2(from_pcs$embedding, curve[c("u", "v")]), 0.995)
})

test_that("invalid arguments are errors", {
  expect_snapshot(error = TRUE, {
    manifold_sculpting(x, ndim = 3)
    manifold_sculpting(x, ndim = 0)
    manifold_sculpting(x, neighbors = 500)
    manifold_sculpting(x, squishing_rate = 1)
    manifold_sculpting(x, nstart = 0)
    manifold_sculpting(x, init = "isomap")
    manifold_sculpting(x, init = x[1:10, ])
    manifold_sculpting(x, nn_method = "ball_tree")
    manifold_sculpting(x, verbose = NA)
    manifold_sculpting(letters)
    manifold_sculpting(data.frame(a = 1:3, b = letters[1:3]))
    manifold_sculpting(rbind(x, NA))
    manifold_sculpting(matrix(1, 10, 3), neighbors = 3)
  })
})

test_that("`verbose` reports the progress", {
  set.seed(1)
  expect_output(
    manifold_sculpting(x, neighbors = 14, max_iter = 50, verbose = TRUE),
    "iteration 50"
  )
})

test_that("print() describes the fit", {
  expect_output(print(fit), "Manifold Sculpting embedding")
  expect_output(print(fit), "500 samples, 3 variables to 2 dimensions")
})
