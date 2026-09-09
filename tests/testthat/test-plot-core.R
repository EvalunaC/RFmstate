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
  with_plot_device(expect_silent(plot(
    aj, type = "cumulative_hazard", states = "B",
    xlab = "Study day", ylab = "Cumulative hazard"
  )))
  with_plot_device(expect_warning(
    plot(aj, type = "transition_intensity"), "deprecated"
  ))
  with_plot_device(expect_warning(plot(aj, ci = TRUE), "unavailable"))
  expect_error(plot(aj, ci = NA), "'ci'.*nonmissing logical")
  expect_error(plot(aj, ci = 1), "'ci'.*nonmissing logical")
  expect_error(plot(aj, states = "unknown"), "Unknown state")
  expect_error(plot(aj, states = character()), "nonempty character")
  expect_error(plot(aj, states = c("B", "B")), "unique state")
  expect_error(
    plot(aj, type = "cumulative_hazard", states = "A"),
    "No transitions lead"
  )
  expect_equal(
    RFmstate:::.filter_plot_transitions(ms$trans_list, "B")$to,
    "B"
  )
})

test_that("prediction state filtering and diagnostic plotting contracts work", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  entry <- array(
    c(1, 0, 0.5, 0.5),
    dim = c(1, 1, 2, 2),
    dimnames = list("profile_1", "A", c("A", "B"), c("0", "1"))
  )
  pred <- structure(list(
    time = c(0, 1), P = entry, entry_prob = entry,
    state_occ = array(c(1, 0, 0.5, 0.5), dim = c(1, 2, 2)),
    structure = ms, n_subjects = 1, start_state = "A"
  ), class = "rfmstate_pred")
  with_plot_device(expect_silent(plot(
    pred, states = "B", xlab = "Duration", ylab = "State probability"
  )))
  with_plot_device(expect_silent(plot(
    pred, type = "transition_prob", states = "B"
  )))
  with_plot_device(expect_silent(plot(
    pred, type = "transition_prob", subject = 0, states = "B"
  )))
  expect_error(plot(pred, states = "Z"), "Unknown state")
  expect_error(plot(pred, subject = -1), "subject")
  expect_error(plot(pred, subject = 2), "exceeds")

  diagnostics <- structure(list(
    concordance = data.frame(transition = "A->B", c_index = 0.7),
    brier = NULL
  ), class = "rfmstate_diag")
  with_plot_device(expect_silent(plot(diagnostics, type = "concordance")))
  expect_error(plot(diagnostics, type = "brier"), "No patient-level")
})

test_that("importance and transition-diagram plots expose axis labels", {
  importance_object <- structure(list(
    importance_matrix = matrix(
      c(0.2, -0.1), nrow = 2,
      dimnames = list(c("x", "z"), "A->B")
    )
  ), class = "rfmstate_importance")
  with_plot_device(expect_silent(plot(
    importance_object, type = "barplot",
    xlab = "OOB loss increase", ylab = "Predictor"
  )))
  with_plot_device(expect_silent(plot(
    importance_object, type = "heatmap",
    xlab = "Edge", ylab = "Predictor"
  )))

  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  with_plot_device(expect_silent(plot_transition_diagram(
    ms, xlab = "Graph layer", ylab = "Display position"
  )))
  expect_error(plot_transition_diagram(list()), "mstate_structure")
})
