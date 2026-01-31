suppressPackageStartupMessages({
  library(testthat)
})

source("app.R")

test_that("haplotype frequencies sum to 1", {
  cfg <- list(
    p = c(MF = 0.1, Mf = 0.2, Mz = 0.1, mF = 0.2, mf = 0.2, mz = 0.2),
    s_M = 0.02,
    s_m = 0.4,
    u_m = 0.01,
    u_f = 0.03,
    u_z = 0.02,
    z_reduction = 0.5,
    max_generations = 25,
    tol = 1e-8
  )

  df <- simulate_recursion(cfg)
  sums <- rowSums(df[, c("MF", "Mf", "Mz", "mF", "mf", "mz")])
  expect_true(all(abs(sums - 1) < 1e-10))
})
