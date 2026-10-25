test_that("the Swiss roll follows the paper", {
  set.seed(1)
  roll <- sim_swiss_roll(200)
  expect_s3_class(roll, "tbl_df")
  expect_named(roll, c("x", "y", "z", "u", "v"))
  t <- 8 * (0:199) / 200 + 2
  expect_equal(roll$x, t * sin(t))
  expect_equal(roll$z, t * cos(t))
  expect_identical(roll$v, roll$y)
  expect_true(all(roll$y >= -6 & roll$y < 6))
  # u is the arc length along the roll: its derivative is the speed
  speed <- sqrt(t^2 + 1)
  expect_equal(diff(roll$u) / diff(t), (speed[-1] + speed[-200]) / 2, tolerance = 1e-3)
})

test_that("the S-curve follows the paper", {
  set.seed(1)
  curve <- sim_s_curve(300)
  t <- (2.2 * (0:299) - 0.1) * pi / 300
  expect_equal(curve$x, t)
  expect_equal(curve$y, sin(t))
  expect_identical(curve$v, curve$z)
  arc <- vapply(t, function(s) integrate(function(w) sqrt(cos(w)^2 + 1), 0, s)$value, 0)
  expect_equal(curve$u, arc, tolerance = 1e-6)
})

test_that("noise is added to the coordinates", {
  set.seed(1)
  clean <- sim_swiss_roll(100)
  set.seed(1)
  noisy <- sim_swiss_roll(100, noise = 0.5)
  expect_identical(noisy$u, clean$u)
  expect_identical(noisy$v, clean$v)
  expect_gt(sd(noisy$x - clean$x), 0.3)
  expect_snapshot(error = TRUE, {
    sim_swiss_roll(10, noise = -1)
    sim_s_curve(1.5)
  })
})
