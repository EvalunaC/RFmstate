test_that("define_multistate creates valid structure", {
  ms <- define_multistate(
    state_names = c("A", "B", "C"),
    absorbing = "C",
    transitions = list(A = c("B", "C"), B = c("C"))
  )

  expect_s3_class(ms, "mstate_structure")
  expect_equal(ms$n_states, 3)
  expect_equal(ms$n_transitions, 3)
  expect_equal(ms$absorbing, "C")
  expect_equal(ms$transient, c("A", "B"))

  # Transition matrix
  expect_equal(ms$trans_matrix["A", "B"], 1L)
  expect_equal(ms$trans_matrix["A", "C"], 2L)
  expect_equal(ms$trans_matrix["B", "C"], 3L)
  expect_true(is.na(ms$trans_matrix["C", "A"]))
  expect_true(is.na(ms$trans_matrix["A", "A"]))
})

test_that("clinical_states creates standard structure", {
  ms <- clinical_states()
  expect_s3_class(ms, "mstate_structure")
  expect_equal(ms$n_states, 6)
  expect_equal(ms$absorbing, "Death")
  expect_equal(length(ms$transient), 5)
  expect_true(ms$n_transitions >= 10)
})

test_that("define_multistate validates inputs", {
  expect_error(
    define_multistate(c("A"), "A", list()),
    "at least 2 states"
  )
  expect_error(
    define_multistate(c("A", "B"), "X", list(A = "B")),
    "absorbing states must be in"
  )
  expect_error(
    define_multistate(c("A", "B"), "B", list(A = c("B", "X"))),
    "Unknown destination"
  )
  expect_error(
    define_multistate(c("A", "B", "C"), "C", list(A = c("B", "C"))),
    "missing from 'transitions'"
  )
})

test_that("print.mstate_structure works", {
  ms <- clinical_states()
  expect_output(print(ms), "Multistate Structure")
  expect_output(print(ms), "Transitions:")
})
