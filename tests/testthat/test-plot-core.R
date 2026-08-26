with_plot_device <- function(code) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  force(code)
}

test_that("AJ plots accept axis labels and hazard-increment terminology", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  dat <- data.frame(id = 1:6, x = 1:6, time_B = 1:6, censor = NA_real_)
  long <- prepare_data(dat, "id", ms, list(B = "time_B"), "censor", "x")
  aj <- aalen_johansen(long)
  with_plot_device(expect_silent(plot(
    aj, type = "state_occupation", xlab = "Study day", ylab = "Probability"
  )))
  with_plot_device(expect_silent(plot(aj, type = "hazard_increment")))
  with_plot_device(expect_warning(
    plot(aj, type = "transition_intensity"), "deprecated"
  ))
  expect_error(plot(aj, states = "unknown"), "Unknown state")
})

test_that("prediction state filtering and diagnostic plotting contracts work", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  entry <- array(0, dim = c(1, 2, 2, 2))
  entry[1, , , 1] <- diag(2)
  entry[1, 1, , 2] <- c(0.5, 0.5)
  entry[1, 2, , 2] <- c(0, 1)
  pred <- structure(list(
    time = c(0, 1), P = entry, entry_prob = entry,
    state_occ = array(c(1, 0, 0.5, 0.5), dim = c(1, 2, 2)),
    structure = ms, n_subjects = 1, start_state = "A"
  ), class = "rfmstate_pred")
  with_plot_device(expect_silent(plot(
    pred, states = "B", xlab = "Duration", ylab = "State probability"
  )))
  expect_error(plot(pred, states = "Z"), "Unknown state")

  diagnostics <- structure(list(
    concordance = data.frame(transition = "A->B", c_index = 0.7),
    brier = NULL
  ), class = "rfmstate_diag")
  with_plot_device(expect_silent(plot(diagnostics, type = "concordance")))
  expect_error(plot(diagnostics, type = "brier"), "No patient-level")
})
