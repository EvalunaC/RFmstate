matrix_exponential_eigen <- function(matrix, time) {
  decomposition <- eigen(matrix)
  value <- decomposition$vectors %*%
    diag(exp(decomposition$values * time), nrow = nrow(matrix)) %*%
    solve(decomposition$vectors)
  Re(value)
}

test_that("branching and converging exponential DAG matches matrix exponential", {
  ms <- define_multistate(
    c("A", "B", "C", "D"), "D",
    list(A = c("B", "C"), B = "D", C = "D")
  )
  rates <- c("A->B" = 0.15, "A->C" = 0.10,
             "B->D" = 0.40, "C->D" = 0.30)
  grid <- seq(0, 3, length.out = 6001)
  curves <- lapply(rates, function(rate) {
    data.frame(time = grid, hazard = rate * grid)
  })
  times <- c(0, 0.5, 1, 2, 3)
  result <- compute_trans_prob(
    curves, ms, times = times, target_grid_points = 1024,
    max_grid_points = 16385, grid_tol = 1e-4
  )
  generator <- matrix(
    c(-0.25, 0.15, 0.10, 0,
       0, -0.40, 0, 0.40,
       0, 0, -0.30, 0.30,
       0, 0, 0, 0),
    nrow = 4, byrow = TRUE
  )
  truth <- t(vapply(times, function(time) {
    matrix_exponential_eigen(generator, time)[1, ]
  }, numeric(4)))
  colnames(truth) <- ms$state_names
  expect_lt(max(abs(result$state_occ - truth)), 1e-3)
})

test_that("non-exponential duration hazards match independent quadrature", {
  ms <- define_multistate(c("A", "B", "C"), "C", list(A = "B", B = "C"))
  shape_a <- 1.4
  scale_a <- 2
  shape_b <- 1.7
  scale_b <- 1.5
  grid <- seq(0, 3, length.out = 6001)
  curves <- list(
    "A->B" = data.frame(
      time = grid, hazard = (grid / scale_a)^shape_a
    ),
    "B->C" = data.frame(
      time = grid, hazard = (grid / scale_b)^shape_b
    )
  )
  times <- c(0, 0.5, 1, 2, 3)
  result <- compute_trans_prob(
    curves, ms, times = times, target_grid_points = 2048,
    max_grid_points = 32769, grid_tol = 1e-4
  )
  truth <- t(vapply(times, function(time) {
    survival_a <- exp(-(time / scale_a)^shape_a)
    occupied_b <- if (time == 0) 0 else stats::integrate(
      function(entry) {
        density_a <- shape_a / scale_a *
          (entry / scale_a)^(shape_a - 1) *
          exp(-(entry / scale_a)^shape_a)
        survival_b <- exp(-((time - entry) / scale_b)^shape_b)
        density_a * survival_b
      },
      lower = 0, upper = time,
      rel.tol = 1e-11, abs.tol = 1e-12
    )$value
    c(A = survival_a, B = occupied_b, C = 1 - survival_a - occupied_b)
  }, numeric(3)))
  expect_lt(max(abs(result$state_occ - truth)), 1e-3)
})

test_that("controlled zero-hazard edge remains zero without masking an edge", {
  ms <- define_multistate(
    c("A", "B", "C"), c("B", "C"), list(A = c("B", "C"))
  )
  grid <- seq(0, 2, length.out = 2001)
  result <- compute_trans_prob(
    list(
      "A->B" = data.frame(time = grid, hazard = 0.2 * grid),
      "A->C" = data.frame(time = grid, hazard = 0 * grid)
    ),
    ms, times = c(0, 1, 2), target_grid_points = 512,
    max_grid_points = 8193, grid_tol = 1e-5
  )
  expect_equal(unname(result$state_occ[, "C"]), c(0, 0, 0), tolerance = 1e-12)
  expect_equal(unname(result$state_occ[, "B"]), 1 - exp(-0.2 * c(0, 1, 2)),
               tolerance = 1e-3)
})

test_that("competing-absorbing AJ agrees with survival multistate estimator", {
  ms <- define_multistate(
    c("A", "B", "C"), c("B", "C"), list(A = c("B", "C"))
  )
  dat <- data.frame(
    id = 1:6,
    x = 1:6,
    time_B = c(1, NA, 3, NA, NA, NA),
    time_C = c(NA, 2, NA, NA, 5, NA),
    censor = c(NA, NA, NA, 4, NA, 6)
  )
  long <- prepare_data(
    dat, "id", ms, list(B = "time_B", C = "time_C"), "censor", "x"
  )
  aj <- aalen_johansen(long)
  event <- factor(
    c("B", "C", "B", "censor", "C", "censor"),
    levels = c("censor", "B", "C")
  )
  reference <- survival::survfit(
    survival::Surv(1:6, event) ~ 1,
    data = data.frame(id = 1:6), id = id
  )
  reference_rows <- match(aj$time, reference$time)
  reference_probability <- reference$pstate[reference_rows, , drop = FALSE]
  colnames(reference_probability) <- c("A", "B", "C")
  expect_equal(aj$state_occ, reference_probability, tolerance = 1e-12)
})

test_that("censored full-state IPCW contribution matches direct calculation", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  training <- data.frame(
    id = 1:4, from = factor("A", levels = c("A", "B")),
    to = factor(c("B", NA, "B", NA), levels = c("A", "B")),
    Tstart = 0, Tstop = c(1, 2, 3, 4),
    status = c(1L, 0L, 1L, 0L), trans_id = c(1L, NA, 1L, NA),
    duration = c(1, 2, 3, 4)
  )
  class(training) <- c("msdata", "data.frame")
  attr(training, "structure") <- ms
  km <- RFmstate:::.fit_censoring_km(training)

  validation <- data.frame(
    id = 1:3, from = factor("A", levels = c("A", "B")),
    to = factor(c("B", NA, "B"), levels = c("A", "B")),
    Tstart = 0, Tstop = c(1, 1.5, 3),
    status = c(1L, 0L, 1L), trans_id = c(1L, NA, 1L),
    duration = c(1, 1.5, 3)
  )
  prediction <- list(
    structure = ms,
    initial_state = "A",
    state_occ = array(
      c(0.1, 0.5, 0.7, 0.9, 0.5, 0.3),
      dim = c(3, 2, 1),
      dimnames = list(NULL, c("A", "B"), "2")
    )
  )
  scored <- RFmstate:::.score_full_state_fold(
    validation, prediction, km, eval_times = 2, g_min = 0.05
  )
  direct <- 1 * (0.1^2 + (1 - 0.9)^2) +
    1.5 * ((1 - 0.7)^2 + 0.3^2)
  expect_equal(scored$contribution_sum, direct, tolerance = 1e-12)
  expect_equal(scored$n_evaluable, 2L)
  expect_equal(direct / 3, 0.0966666666666667, tolerance = 1e-12)
})

test_that("invalid OOB backend results fail explicitly", {
  expect_error(
    RFmstate:::.validate_oob_result(
      list(prediction.error = NaN, inbag.counts = list(c(0, 1))), 2, "A->B"
    ),
    "finite valid OOB"
  )
  expect_error(
    RFmstate:::.validate_oob_result(
      list(prediction.error = 0.3, inbag.counts = list(c(1, 1))), 2, "A->B"
    ),
    "no verified OOB"
  )
})
