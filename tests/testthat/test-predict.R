# The points of sim_s_curve() are ordered along the curve, so the split is
# random
data <- local({
  set.seed(20)
  curve <- sim_s_curve(700)
  curve[sample(700), ]
})
train <- data[1:500, ]
test <- data[501:700, ]
fit <- local({
  set.seed(1)
  manifold_sculpting(train[c("x", "y", "z")], neighbors = 14)
})

# Maps an embedding onto the intrinsic coordinates, with the Procrustes
# transformation of the training embedding
to_truth <- function(fit, truth) {
  y <- fit$embedding
  my <- colMeans(y)
  mt <- colMeans(truth)
  s <- svd(crossprod(sweep(y, 2, my), sweep(as.matrix(truth), 2, mt)))
  b <- sum(s$d) / sum(sweep(y, 2, my)^2)
  function(z) sweep(b * sweep(z, 2, my) %*% s$u %*% t(s$v), 2, mt, "+")
}

test_that("the training data are mapped to the embedding", {
  expect_identical(predict(fit, train), fit$embedding)
  expect_identical(predict(fit, train, refine = FALSE), fit$embedding)
})

test_that("new samples are embedded as accurately as the training samples", {
  map <- to_truth(fit, train[c("u", "v")])
  rmse <- function(y, truth) sqrt(mean(rowSums((map(y) - as.matrix(truth))^2)))
  train_rmse <- rmse(fit$embedding, train[c("u", "v")])

  pred <- predict(fit, test)
  expect_identical(dim(pred), c(200L, 2L))
  expect_identical(colnames(pred), c("MS1", "MS2"))
  expect_lt(rmse(pred, test[c("u", "v")]), 2 * train_rmse + 0.01)

  lle <- predict(fit, test, refine = FALSE)
  expect_lt(rmse(lle, test[c("u", "v")]), 2 * train_rmse + 0.01)
})

test_that("refinement corrects the extrapolation of the reconstruction weights", {
  # Train on the first 85% of the curve and embed the next points
  set.seed(21)
  curve <- sim_s_curve(600)
  first <- curve[curve$u < quantile(curve$u, 0.85), ]
  beyond <- curve[curve$u >= quantile(curve$u, 0.85), ][1:30, ]
  set.seed(1)
  ordered_fit <- manifold_sculpting(first[c("x", "y", "z")], neighbors = 14)
  map <- to_truth(ordered_fit, first[c("u", "v")])
  rmse <- function(y) sqrt(mean(rowSums((map(y) - as.matrix(beyond[c("u", "v")]))^2)))
  refined <- rmse(predict(ordered_fit, beyond))
  lle <- rmse(predict(ordered_fit, beyond, refine = FALSE))
  expect_lt(refined, lle)
  expect_lt(refined, 0.1)
})

test_that("columns are matched by name", {
  shuffled <- test[c("z", "u", "x", "y")]
  expect_identical(predict(fit, shuffled), predict(fit, test))
  expect_identical(predict(fit, as.matrix(shuffled)), predict(fit, test))
  expect_snapshot(error = TRUE, predict(fit, test[c("x", "y")]))
})

test_that("columns are matched by position without names", {
  set.seed(1)
  unnamed <- manifold_sculpting(unname(as.matrix(train[c("x", "y", "z")])), neighbors = 14)
  new <- unname(as.matrix(test[c("x", "y", "z")]))
  expect_identical(unname(predict(unnamed, new)), unname(predict(fit, test)))
  expect_snapshot(error = TRUE, predict(unnamed, new[, 1:2]))
})

test_that("edge cases", {
  expect_identical(dim(predict(fit, test[0, ])), c(0L, 2L))
  expect_snapshot(error = TRUE, {
    predict(fit, test, refine = NA)
    predict(fit, test, foo = 1)
  })
})
