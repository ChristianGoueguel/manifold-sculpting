# R^2 of the intrinsic coordinates `truth` predicted from the embedding `y`
# by the best rotation, scaling and translation (as in the paper).
procrustes_r2 <- function(y, truth) {
  y <- scale(as.matrix(y), scale = FALSE)
  truth <- scale(as.matrix(truth), scale = FALSE)
  s <- svd(crossprod(y, truth))
  b <- sum(s$d) / sum(y^2)
  fitted <- b * y %*% s$u %*% t(s$v)
  1 - sum((truth - fitted)^2) / sum(truth^2)
}

# Brute-force neighbors, computed in R, with ties broken by the index
r_knn <- function(ref, query, k, exclude_self) {
  d2 <- 0
  for (j in seq_len(ncol(ref))) {
    d2 <- d2 + outer(query[, j], ref[, j], "-")^2
  }
  if (exclude_self) diag(d2) <- Inf
  index <- t(apply(d2, 1, function(r) order(r, seq_along(r))[seq_len(k)]))
  list(index = index, distance = sqrt(matrix(d2[cbind(rep(seq_len(nrow(query)), k), c(index))], ncol = k)))
}
