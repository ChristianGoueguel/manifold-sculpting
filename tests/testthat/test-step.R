data <- local({
  set.seed(30)
  curve <- sim_s_curve(400)
  curve$label <- factor(rep(c("a", "b"), 200))
  curve
})
train <- data[1:300, ]
test <- data[301:400, ]

rec <- recipes::recipe(label ~ x + y + z, data = train)

test_that("the step adds the embedding and removes the predictors", {
  trained <- rec |>
    step_manifold_sculpting(x, y, z, neighbors = 14, seed = 1) |>
    prep()

  res <- trained$steps[[1]]$res
  expect_s3_class(res, "manifold_sculpting")

  # The training set gets the embedding of the algorithm
  juiced <- bake(trained, new_data = NULL)
  expect_named(juiced, c("label", "MS1", "MS2"))
  expect_equal(as.matrix(juiced[c("MS1", "MS2")]), unname(res$embedding),
               ignore_attr = TRUE)
  expect_equal(bake(trained, new_data = train), juiced)

  # New data are embedded by predict()
  baked <- bake(trained, new_data = test)
  expect_named(baked, c("label", "MS1", "MS2"))
  expect_equal(as.matrix(baked[c("MS1", "MS2")]),
               predict(res, test), ignore_attr = TRUE)
})

test_that("the step matches manifold_sculpting() with the same seed", {
  trained <- rec |>
    step_manifold_sculpting(recipes::all_numeric_predictors(),
                            num_comp = 2, neighbors = 14, seed = 7) |>
    prep()
  direct <- withr::with_seed(7, manifold_sculpting(train[c("x", "y", "z")], neighbors = 14))
  expect_identical(trained$steps[[1]]$res$embedding, direct$embedding)
})

test_that("arguments are passed on", {
  trained <- rec |>
    step_manifold_sculpting(x, y, z, num_comp = 1, neighbors = 10,
                            squishing_rate = 0.95, prefix = "dim",
                            options = list(max_iter = 30, nstart = 2), seed = 2) |>
    prep()
  res <- trained$steps[[1]]$res
  expect_identical(res$ndim, 1L)
  expect_identical(res$neighbors, 10L)
  expect_identical(res$squishing_rate, 0.95)
  expect_length(res$final_error, 2)
  expect_identical(res$iterations, 30L)
  expect_named(bake(trained, new_data = test), c("label", "dim1"))
})

test_that("keep_original_cols keeps the predictors", {
  trained <- rec |>
    step_manifold_sculpting(x, y, z, neighbors = 14, keep_original_cols = TRUE, seed = 1) |>
    prep()
  expect_named(bake(trained, new_data = test), c("x", "y", "z", "label", "MS1", "MS2"))
})

test_that("num_comp is reduced below the number of columns", {
  expect_snapshot(
    trained <- rec |>
      step_manifold_sculpting(x, y, neighbors = 14, num_comp = 3, seed = 1,
                              options = list(max_iter = 20)) |>
      prep()
  )
  expect_equal(trained$steps[[1]]$num_comp, 1)
  expect_named(bake(trained, new_data = test), c("z", "label", "MS1"))
})

test_that("num_comp = 0 and empty selections leave the data unchanged", {
  trained <- rec |>
    step_manifold_sculpting(x, y, z, num_comp = 0) |>
    prep()
  expect_null(trained$steps[[1]]$res)
  expect_identical(bake(trained, new_data = test), test[c("x", "y", "z", "label")])

  trained <- rec |>
    step_manifold_sculpting() |>
    prep()
  expect_identical(bake(trained, new_data = test), test[c("x", "y", "z", "label")])
})

test_that("baking zero rows works", {
  trained <- rec |>
    step_manifold_sculpting(x, y, z, neighbors = 14, seed = 1, options = list(max_iter = 20)) |>
    prep()
  baked <- bake(trained, new_data = test[0, ])
  expect_named(baked, c("label", "MS1", "MS2"))
  expect_identical(nrow(baked), 0L)
})

test_that("missing columns at bake time are an error", {
  trained <- rec |>
    step_manifold_sculpting(x, y, z, neighbors = 14, seed = 1, options = list(max_iter = 20)) |>
    prep()
  expect_snapshot(error = TRUE, bake(trained, new_data = test[c("x", "y", "label")]))
})

test_that("invalid arguments are errors", {
  expect_snapshot(error = TRUE, {
    rec |> step_manifold_sculpting(x, y, z, options = list(foo = 1)) |> prep()
    rec |> step_manifold_sculpting(x, y, z, squishing_rate = 2) |> prep()
    rec |> step_manifold_sculpting(x, y, z, neighbors = -1) |> prep()
    rec |> step_manifold_sculpting(x, y, z, prefix = 1) |> prep()
    rec |> step_manifold_sculpting(label) |> prep()
  })
})

test_that("printing", {
  step <- rec |> step_manifold_sculpting(x, y, z, neighbors = 14, seed = 1,
                                         options = list(max_iter = 20))
  expect_snapshot(print(step))
  expect_snapshot(prep(step))
  expect_snapshot(rec |> step_manifold_sculpting(x, y, z, num_comp = 0))
})

test_that("tidy() lists the selected columns", {
  step <- rec |> step_manifold_sculpting(x, y, z, neighbors = 14, seed = 1,
                                         options = list(max_iter = 20), id = "ms")
  untrained <- tidy(step, number = 1)
  expect_identical(untrained$terms, c("x", "y", "z"))
  expect_identical(untrained$id, rep("ms", 3))

  trained <- tidy(prep(step), number = 1)
  expect_identical(trained$terms, c("x", "y", "z"))
  expect_identical(trained$id, rep("ms", 3))
})

test_that("tunable parameters", {
  step <- rec |> step_manifold_sculpting(x, y, z)
  params <- tunable(step$steps[[1]])
  expect_identical(params$name, c("num_comp", "neighbors", "squishing_rate"))
  expect_identical(params$component, rep("step_manifold_sculpting", 3))

  rate <- squishing_rate()
  expect_s3_class(rate, "quant_param")
  expect_identical(rate$range, list(lower = 0.95, upper = 0.995))
})

test_that("required packages", {
  step <- rec |> step_manifold_sculpting(x, y, z)
  expect_identical(recipes::required_pkgs(step), c("recipes", "manifoldsculpting"))
})

test_that("the step works in a workflow", {
  skip_if_not_installed("workflows")
  skip_if_not_installed("parsnip")
  wf <- workflows::workflow() |>
    workflows::add_recipe(
      recipes::recipe(z ~ x + y + u, data = train) |>
        step_manifold_sculpting(x, y, u, num_comp = 2, neighbors = 14, seed = 1,
                                options = list(max_iter = 50))
    ) |>
    workflows::add_model(parsnip::linear_reg())
  fitted <- parsnip::fit(wf, data = train)
  pred <- predict(fitted, new_data = test)
  expect_identical(nrow(pred), 100L)
})
