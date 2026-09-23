context("N_alloc_Power function")

## Regression tests: combn() expands a length-1 index vector to seq_len(n),
## which previously produced bogus allocations (repeated/missing cluster
## sizes), and the per-allocation powers were kept as a list, breaking
## quantile().

test_that("exhaustive enumeration yields valid permutations only", {
  sizes <- c(12, 8, 10, 9, 14)
  ap <- N_alloc_Power(N_pool = sizes, Cl = rep(1, 5),
                      mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33,
                      keep = "all", verbose = FALSE)
  expect_equal(ap$summary[["n_allocations"]], 120)  # 5! index-partitions
  expect_equal(ap$summary[["n_unique"]], 120)
  expect_equal(sum(ap$weights), 120)
  expect_true(all(apply(ap$N_allocations, 1,
                        function(r) identical(sort(r), sort(sizes)))))
  expect_identical(sort(ap$best), sort(sizes))
  expect_identical(sort(ap$worst), sort(sizes))
})

test_that("powers are atomic and consistent with direct glsPower calls", {
  sizes <- c(12, 8, 10, 9, 14)
  ap <- N_alloc_Power(N_pool = sizes, Cl = rep(1, 5),
                      mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33,
                      verbose = FALSE)
  expect_is(ap$power, "numeric")
  p_best <- glsPower(Cl = rep(1, 5), N = ap$best,
                     mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33, verbose = 0)
  p_worst <- glsPower(Cl = rep(1, 5), N = ap$worst,
                      mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33, verbose = 0)
  expect_equal(max(ap$power), p_best)
  expect_equal(min(ap$power), p_worst)
})

test_that("multi-cluster sequences and duplicate sizes are handled", {
  ## 8 clusters, 4 sequences of 2: 8!/(2!^4) = 2520 index-partitions
  sizes <- c(20, 30, 25, 15, 40, 10, 35, 22)
  ap <- N_alloc_Power(N_pool = sizes, Cl = c(2, 2, 2, 2),
                      mu0 = 0.2, mu1 = 0.35, tau = 0.5,
                      family = "binomial", verbose = FALSE)
  expect_equal(ap$summary[["n_allocations"]], 2520)

  ## duplicate sizes are de-duplicated and weighted
  apd <- N_alloc_Power(N_pool = c(10, 10, 20, 20), Cl = c(2, 2),
                       mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33,
                       verbose = FALSE)
  expect_equal(apd$summary[["n_allocations"]], 6)
  expect_equal(apd$summary[["n_unique"]], 3)
  expect_equal(sum(apd$weights), 6)
})

test_that("Monte Carlo mode works", {
  set.seed(123)
  ap <- N_alloc_Power(N_pool = c(12, 8, 10, 9, 14), Cl = rep(1, 5),
                      mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33,
                      n_MC = 20, verbose = FALSE)
  expect_is(ap$power, "numeric")
  expect_length(ap$power, 20)
  expect_equal(ap$weights, rep(1, 20))
  expect_equal(ap$summary[["n_allocations"]], 20)
})
