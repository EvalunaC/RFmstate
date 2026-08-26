simple_structure <- function() {
  define_multistate(c("A", "B", "C"), "C", list(A = "B", B = "C"))
}

simple_wide <- function() {
  data.frame(
    id = 1:4,
    x = c(0, 1, 0, 1),
    time_B = c(1, 2, NA, 1.5),
    time_C = c(3, 4, NA, 3.5),
    censor = c(NA, NA, 5, NA)
  )
}

prepare_simple <- function(data = simple_wide()) {
  prepare_data(data, "id", simple_structure(),
               list(B = "time_B", C = "time_C"),
               "censor", "x")
}

test_that("graph validation enforces a single-root DAG and preserves display order", {
  display_backward <- define_multistate(
    c("C", "A", "B"), "C", list(A = "B", B = "C")
  )
  expect_equal(display_backward$state_names, c("C", "A", "B"))
  expect_equal(display_backward$topological_order, c("A", "B", "C"))
  expect_equal(display_backward$initial_state, "A")

  expect_error(
    define_multistate(c("A", "B"), "B", list(A = c("B", "B"))),
    "Duplicated directed edge"
  )
  expect_error(
    define_multistate(
      c("A", "B", "C", "D"), "D",
      list(A = "B", B = "C", C = c("B", "D"))
    ),
    "cycles"
  )
  expect_error(
    define_multistate(c("A", "B", "C"), "C", list(A = "C", B = "C")),
    "exactly one common initial"
  )
})

test_that("prepare_data stores contract metadata and head returns rows", {
  long <- prepare_simple()
  expect_equal(attr(long, "initial_state"), "A")
  expect_equal(attr(long, "covariates"), "x")
  expect_equal(attr(long, "metadata")$n_subjects, 4)
  expect_equal(attr(long, "metadata")$display_order, c("A", "B", "C"))
  expect_s3_class(head(long), "data.frame")
  expect_false(inherits(head(long), "msdata"))
  expect_true(all(long$duration > 0))
})

test_that("prepare_data rejects invalid IDs, maps, times, and paths", {
  dat <- simple_wide()
  dat$id[2] <- dat$id[1]
  expect_error(prepare_simple(dat), "IDs must be unique")

  dat <- simple_wide()
  dat$x[1] <- NA
  expect_error(prepare_simple(dat), "contains missing")

  expect_error(
    prepare_data(simple_wide(), "id", simple_structure(),
                 list(B = "time_B"), "censor", "x"),
    "every noninitial state"
  )

  dat <- simple_wide()
  dat$time_C[1] <- dat$time_B[1]
  expect_error(prepare_simple(dat), "simultaneous")

  dat <- simple_wide()
  dat$time_B[1] <- 0
  expect_error(prepare_simple(dat), "nonincreasing")

  dat <- simple_wide()
  dat$censor[1] <- 2
  expect_error(prepare_simple(dat), "after censoring")

  dat <- simple_wide()
  dat$time_B[1] <- NA
  expect_error(prepare_simple(dat), "forbidden observed transition")
})

test_that("events after absorption are rejected", {
  ms <- define_multistate(
    c("A", "B", "C"), "C", list(A = c("B", "C"), B = "C")
  )
  dat <- data.frame(id = 1, x = 1, time_B = 3, time_C = 2, censor = NA)
  expect_error(
    prepare_data(dat, "id", ms, list(B = "time_B", C = "time_C"),
                 "censor", "x"),
    "after absorption"
  )
})
