knn <- function(...) manifoldsculpting:::knn_cpp(...)

test_that("the k-d tree and brute-force searches find the exact neighbors", {
  set.seed(1)
  x <- matrix(rnorm(300 * 4), ncol = 4)
  q <- matrix(rnorm(40 * 4), ncol = 4)
  for (use_tree in c(TRUE, FALSE)) {
    res <- knn(x, x, 7, TRUE, use_tree)
    ref <- r_knn(x, x, 7, exclude_self = TRUE)
    expect_identical(res$index, ref$index)
    expect_equal(res$distance, ref$distance, tolerance = 1e-12)

    res <- knn(x, q, 5, FALSE, use_tree)
    ref <- r_knn(x, q, 5, exclude_self = FALSE)
    expect_identical(res$index, ref$index)
    expect_equal(res$distance, ref$distance, tolerance = 1e-12)
  }
})

test_that("ties are broken by the index", {
  # A grid has many equidistant neighbors
  x <- as.matrix(expand.grid(1:12, 1:12, 1:3))
  storage.mode(x) <- "double"
  ref <- r_knn(x, x, 9, exclude_self = TRUE)
  expect_identical(knn(x, x, 9, TRUE, TRUE)$index, ref$index)
  expect_identical(knn(x, x, 9, TRUE, FALSE)$index, ref$index)
})

test_that("duplicated points and high dimensions are handled", {
  set.seed(2)
  x <- matrix(rnorm(100 * 30), ncol = 30)
  x <- rbind(x, x[1:10, ])
  tree <- knn(x, x, 4, TRUE, TRUE)
  brute <- knn(x, x, 4, TRUE, FALSE)
  expect_identical(tree$index, brute$index)
  expect_equal(tree$distance[, 1][1:10], rep(0, 10))
  expect_equal(tree$index[1:10, 1], 101:110)
})

test_that("invalid k is an error", {
  x <- matrix(rnorm(20), ncol = 2)
  expect_error(knn(x, x, 10, TRUE, TRUE), "between 1 and")
  expect_error(knn(x, x, 0, FALSE, TRUE), "between 1 and")
  expect_error(knn(x, x[, 1, drop = FALSE], 2, FALSE, TRUE), "same number of columns")
})

test_that("`nn_method` chooses the search", {
  expect_true(manifoldsculpting:::use_kd_tree("auto", 3))
  expect_false(manifoldsculpting:::use_kd_tree("auto", 50))
  expect_true(manifoldsculpting:::use_kd_tree("kd_tree", 50))
  expect_false(manifoldsculpting:::use_kd_tree("brute_force", 3))
})
