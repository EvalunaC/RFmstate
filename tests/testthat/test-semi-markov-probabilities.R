hazard_curve <- function(rate, horizon = 3, n = 6001L) {
  time <- seq(0, horizon, length.out = n)
  data.frame(time = time, hazard = rate * time)
}

test_that("two-state exponential probabilities match analytic truth", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  times <- c(0, 0.5, 1, 2)
  result <- compute_trans_prob(
    list("A->B" = hazard_curve(0.4)), ms, times = times,
    target_grid_points = 512, max_grid_points = 8193, grid_tol = 1e-5
  )
  expect_lt(max(abs(result$state_occ[, "A"] - exp(-0.4 * times))), 1e-3)
  expect_lt(max(abs(result$state_occ[, "B"] -
                      (1 - exp(-0.4 * times)))), 1e-3)
  expect_equal(unname(result$entry_prob["B", "B", ]), rep(1, length(times)),
               tolerance = 1e-10)
})

test_that("progressive exponential DAG matches analytic truth", {
  ms <- define_multistate(
    c("A", "B", "C"), "C", list(A = c("B", "C"), B = "C")
  )
  a <- 0.2; b <- 0.1; c_rate <- 0.4
  times <- c(0, 0.5, 1, 2)
  result <- compute_trans_prob(
    list(
      "A->B" = hazard_curve(a),
      "A->C" = hazard_curve(b),
      "B->C" = hazard_curve(c_rate)
    ), ms, times = times, target_grid_points = 512,
    max_grid_points = 8193, grid_tol = 1e-5
  )
  p_a <- exp(-(a + b) * times)
  p_b <- a * exp(-c_rate * times) *
    (1 - exp(-(a + b - c_rate) * times)) / (a + b - c_rate)
  p_b[1L] <- 0
  truth <- cbind(A = p_a, B = p_b, C = 1 - p_a - p_b)
  expect_lt(max(abs(result$state_occ - truth)), 1e-3)
})

test_that("competing absorbing destinations and probability invariants hold", {
  ms <- define_multistate(
    c("A", "B", "C"), c("B", "C"), list(A = c("B", "C"))
  )
  a <- 0.3; b <- 0.2
  times <- c(0, 0.5, 1, 2)
  result <- compute_trans_prob(
    list("A->B" = hazard_curve(a), "A->C" = hazard_curve(b)),
    ms, times = times, target_grid_points = 512,
    max_grid_points = 8193, grid_tol = 1e-5
  )
  exit <- 1 - exp(-(a + b) * times)
  truth <- cbind(A = 1 - exit, B = a / (a + b) * exit,
                 C = b / (a + b) * exit)
  expect_lt(max(abs(result$state_occ - truth)), 1e-3)
  expect_true(all(is.finite(result$entry_prob)))
  expect_gte(min(result$entry_prob), -1e-10)
  expect_lte(max(result$entry_prob), 1 + 1e-10)
  expect_lt(max(abs(apply(result$entry_prob, c(1, 3), sum) - 1)), 1e-8)
  expect_equal(unname(result$entry_prob[, , 1]), diag(3), tolerance = 1e-10)
})

test_that("solver validates edges, hazards, support, and output-time spacing", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  curve <- hazard_curve(0.4, horizon = 2)
  expect_error(compute_trans_prob(list(), ms, times = 1), "uniquely named")
  bad <- curve
  bad$hazard[100] <- bad$hazard[99] - 0.01
  expect_error(compute_trans_prob(list("A->B" = bad), ms, times = 1),
               "decreases")
  expect_error(compute_trans_prob(list("A->B" = curve), ms, times = 3),
               "exceeds conservative support")
  common <- c(0, 0.5, 1, 1.5, 2)
  p1 <- compute_trans_prob(
    list("A->B" = curve), ms, times = common,
    target_grid_points = 256, max_grid_points = 4097
  )
  p2 <- compute_trans_prob(
    list("A->B" = curve), ms, times = sort(unique(c(common, seq(0, 2, .1)))),
    target_grid_points = 256, max_grid_points = 4097
  )
  index <- match(common, p2$time)
  expect_equal(p1$state_occ, p2$state_occ[index, ], tolerance = 1e-8)
  expect_error(compute_trans_prob(list("A->B" = curve), ms, s = 1,
                                  times = 1), "Only s = 0")
})

test_that("FFT solver agrees with retained direct convolution reference", {
  ms <- define_multistate(
    c("A", "B", "C"), "C", list(A = c("B", "C"), B = "C")
  )
  curves <- list(
    "A->B" = hazard_curve(0.2, 1, 1001),
    "A->C" = hazard_curve(0.1, 1, 1001),
    "B->C" = hazard_curve(0.3, 1, 1001)
  )
  grid <- seq(0, 1, length.out = 129)
  fft <- RFmstate:::.semi_markov_once(curves, ms, grid, engine = "fft")$prob
  direct <- RFmstate:::.semi_markov_once(curves, ms, grid, engine = "direct")$prob
  expect_equal(fft, direct, tolerance = 1e-10)
})
