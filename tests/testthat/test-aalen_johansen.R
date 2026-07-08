test_that("aalen_johansen computes valid estimates", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 200, structure = ms)
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Responded = "time_Responded",
      Unresponded = "time_Unresponded",
      Stabilized = "time_Stabilized",
      Progressed = "time_Progressed",
      Death = "time_Death"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )

  aj <- aalen_johansen(msdata)

  expect_s3_class(aj, "aj_estimate")
  expect_true(length(aj$time) > 0)

  # State occupation probabilities should sum to ~1
  final_occ <- aj$state_occ[nrow(aj$state_occ), ]
  expect_true(abs(sum(final_occ) - 1) < 0.01)

  # All probabilities non-negative
  expect_true(all(aj$state_occ >= -0.01))

  # Transition probability matrices should have rows summing to ~1
  last_P <- aj$trans_prob[[length(aj$trans_prob)]]
  row_sums <- rowSums(last_P)
  expect_true(all(abs(row_sums - 1) < 0.05))
})

test_that("aalen_johansen handles empty data gracefully", {
  ms <- define_multistate(
    state_names = c("A", "B"),
    absorbing = "B",
    transitions = list(A = "B")
  )

  msdata <- data.frame(
    id = character(0), from = character(0), to = character(0),
    Tstart = numeric(0), Tstop = numeric(0),
    status = integer(0), trans_id = integer(0),
    stringsAsFactors = FALSE
  )
  attr(msdata, "structure") <- ms
  class(msdata) <- c("msdata", "data.frame")

  expect_error(aalen_johansen(msdata), "No events")
})
