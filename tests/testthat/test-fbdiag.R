context("fbdiag function")

## fbdiag assumes symmetric blocks of equal size ####

make_sym <- function(N, seed) {
  set.seed(seed)
  m <- matrix(runif(N * N), N, N)
  m + t(m)
}

test_that("fbdiag matches Matrix::bdiag for symmetric blocks", {
  lmat <- list(make_sym(5, 1), make_sym(5, 2), make_sym(5, 3))
  expect_equal(as.matrix(fbdiag(lmat)),
              as.matrix(Matrix::bdiag(lmat)))
})

test_that("fbdiag works with a single block", {
  lmat <- list(make_sym(4, 1))
  expect_equal(as.matrix(fbdiag(lmat)), lmat[[1]])
})

test_that("fbdiag works with 1x1 blocks", {
  lmat <- lapply(1:6, function(k) matrix(k, 1, 1))
  expect_equal(as.matrix(fbdiag(lmat)), diag(1:6))
})

test_that("fbdiag works with unequal block count and size", {
  lmat <- lapply(1:12, function(k) make_sym(4, k))
  expect_equal(as.matrix(fbdiag(lmat)),
              as.matrix(Matrix::bdiag(lmat)))
})

test_that("fbdiag result is symmetric and off-diagonal blocks are zero", {
  lmat <- list(make_sym(3, 1), make_sym(3, 2))
  cm <- as.matrix(fbdiag(lmat))
  expect_true(isSymmetric(cm))
  expect_true(all(cm[1:3, 4:6] == 0))
  expect_true(all(cm[4:6, 1:3] == 0))
})

test_that("fbdiag diagonal blocks equal input blocks", {
  lmat <- list(make_sym(3, 1), make_sym(3, 2))
  cm <- as.matrix(fbdiag(lmat))
  expect_equal(cm[1:3, 1:3], lmat[[1]])
  expect_equal(cm[4:6, 4:6], lmat[[2]])
})

test_that("fbdiag returns a sparse symmetric Matrix", {
  cm <- fbdiag(list(make_sym(3, 1), make_sym(3, 2)))
  expect_s4_class(cm, "dsCMatrix")
  expect_equal(dim(cm), c(6, 6))
})
