#' Simulate Clinical Trial Multistate Data
#'
#' Generates realistic clinical trial data with covariates and multistate
#' event times for testing and demonstration. The structure must be an acyclic,
#' non-recurrent graph with one common initial state.
#'
#' Transition intensities follow Weibull distributions with covariate effects
#' on the scale parameter. For the default \code{\link{clinical_states}()}
#' structure, transition-specific parameters are calibrated to produce
#' realistic clinical trial trajectories. For custom structures, sensible
#' default parameters are used for all transitions.
#'
#' @param n Integer, number of patients to simulate.
#' @param structure An \code{mstate_structure} object. Defaults to
#'   \code{\link{clinical_states}()}.
#' @param max_followup Numeric, maximum follow-up time (for generating
#'   censoring). Default 365.
#' @param seed Optional integer for reproducibility.
#'
#' @return A data frame in wide format with columns:
#'   \describe{
#'     \item{ID}{Patient identifier (1 to n).}
#'     \item{age}{Continuous, simulated from Normal(60, 12).}
#'     \item{sex}{Binary 0/1.}
#'     \item{BMI}{Continuous, simulated from Normal(26, 5).}
#'     \item{treatment}{Binary 0/1 (balanced arms).}
#'     \item{time_\emph{StateName}}{For each non-initial state in the
#'       structure, the time (days) of entry into that state, or \code{NA}
#'       if the state was not visited. Column names follow the pattern
#'       \code{time_<StateName>} (e.g., \code{time_Death}).}
#'     \item{time_censored}{Days until last follow-up (right censoring
#'       time), or \code{NA} if an absorbing state was reached.}
#'   }
#'
#' @examples
#' set.seed(123)
#' dat <- sim_clinical_data(n = 100)
#' head(dat)
#' summary(dat)
#'
#' @export
sim_clinical_data <- function(n = 500, structure = NULL,
                              max_followup = 365, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(structure)) structure <- clinical_states()
  if (!inherits(structure, "mstate_structure")) {
    stop("'structure' must be an mstate_structure object.")
  }
  if (!is.numeric(n) || length(n) != 1L || n < 1L || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }
  if (!is.numeric(max_followup) || length(max_followup) != 1L ||
      !is.finite(max_followup) || max_followup <= 0) {
    stop("'max_followup' must be positive and finite.")
  }

  # Generate covariates
  age <- rnorm(n, mean = 60, sd = 12)
  sex <- rbinom(n, 1, 0.5)
  bmi <- rnorm(n, mean = 26, sd = 5)
  treatment <- rep(0:1, length.out = n)

  # Initialize result data frame with covariates
  dat <- data.frame(
    ID = seq_len(n),
    age = round(age, 1),
    sex = sex,
    BMI = round(bmi, 1),
    treatment = treatment,
    stringsAsFactors = FALSE
  )

  # Dynamically add time columns for all non-initial states
  non_initial_states <- setdiff(structure$state_names, structure$initial_state)
  for (s in non_initial_states) {
    dat[[paste0("time_", s)]] <- NA_real_
  }
  dat$time_censored <- NA_real_

  # External censoring is drawn before the trajectory and competes with exits.
  cens_time <- runif(n, min = max_followup * 0.3, max = max_followup)

  for (i in seq_len(n)) {
    covs <- c(age[i], sex[i], bmi[i], treatment[i])
    stop_time <- min(cens_time[i], max_followup)
    trajectory <- .sim_patient_trajectory(
      covs, structure, max_followup, censor_time = stop_time
    )

    for (nm in names(trajectory)) {
      col_name <- paste0("time_", nm)
      if (col_name %in% names(dat)) {
        dat[[col_name]][i] <- trajectory[[nm]]
      }
    }

    # Censoring is missing only when absorption was observed before it.
    reached_absorbing <- any(
      names(trajectory) %in% structure$absorbing &
        !is.na(unlist(trajectory[names(trajectory) %in% structure$absorbing]))
    )

    if (reached_absorbing) {
      dat$time_censored[i] <- NA_real_
    } else {
      dat$time_censored[i] <- stop_time
    }
  }

  dat
}

#' Simulate a single patient trajectory
#' @param covs Numeric vector: age, sex, bmi, treatment
#' @param structure mstate_structure
#' @param max_followup Maximum time
#' @noRd
.sim_patient_trajectory <- function(covs, structure, max_followup,
                                    censor_time = max_followup) {
  age <- covs[1]
  sex <- covs[2]
  bmi <- covs[3]
  trt <- covs[4]

  result <- list()
  current_state <- structure$initial_state

  current_time <- 0

  while (!(current_state %in% structure$absorbing) &&
         current_time < max_followup) {
    # Get allowed destinations from current state
    dests <- structure$transitions[[current_state]]
    if (is.null(dests) || length(dests) == 0) break

    # Generate competing event times for each destination
    waits <- vapply(dests, function(dest) {
      .sim_transition_time(current_state, dest, covs)
    }, numeric(1))

    # Find the first event
    min_idx <- which.min(waits)
    event_time <- current_time + waits[min_idx]
    next_state <- dests[min_idx]

    if (event_time >= min(max_followup, censor_time)) break

    result[[next_state]] <- event_time
    current_time <- event_time
    current_state <- next_state
  }

  result
}

#' Simulate transition time from a specific transition
#' @noRd
.sim_transition_time <- function(from, to, covs) {
  age <- covs[1]
  sex <- covs[2]
  bmi <- covs[3]
  trt <- covs[4]

  # Base parameters (shape, scale) for each transition type
  params <- .get_transition_params(from, to)
  shape <- params$shape
  base_scale <- params$scale

  # Covariate effects on scale (log-linear)
  log_scale <- log(base_scale) +
    params$beta_age * (age - 60) / 10 +
    params$beta_sex * sex +
    params$beta_bmi * (bmi - 26) / 5 +
    params$beta_trt * trt

  scale <- exp(log_scale)

  stats::rweibull(1L, shape = shape, scale = scale)
}

#' Get transition-specific simulation parameters
#' @noRd
.get_transition_params <- function(from, to) {
  # Default parameters
  params <- list(
    shape = 1.2, scale = 200,
    beta_age = 0.1, beta_sex = 0.05,
    beta_bmi = 0.05, beta_trt = -0.2
  )

  key <- paste(from, to, sep = "_")
  switch(key,
    "Baseline_Responded" = {
      params$shape <- 1.5
      params$scale <- 60
      params$beta_trt <- -0.5  # Treatment helps response
      params$beta_age <- 0.15
    },
    "Baseline_Unresponded" = {
      params$shape <- 1.3
      params$scale <- 80
      params$beta_trt <- 0.3   # Treatment reduces non-response
      params$beta_age <- 0.1
    },
    "Baseline_Death" = {
      params$shape <- 1.0
      params$scale <- 400
      params$beta_age <- 0.3
      params$beta_bmi <- 0.15
    },
    "Responded_Stabilized" = {
      params$shape <- 1.4
      params$scale <- 90
      params$beta_trt <- -0.3
    },
    "Responded_Progressed" = {
      params$shape <- 1.2
      params$scale <- 180
      params$beta_trt <- 0.2
      params$beta_age <- 0.2
    },
    "Responded_Death" = {
      params$shape <- 1.0
      params$scale <- 500
      params$beta_age <- 0.3
    },
    "Unresponded_Stabilized" = {
      params$shape <- 1.3
      params$scale <- 120
      params$beta_trt <- -0.2
    },
    "Unresponded_Progressed" = {
      params$shape <- 1.5
      params$scale <- 100
      params$beta_trt <- 0.3
      params$beta_age <- 0.15
    },
    "Unresponded_Death" = {
      params$shape <- 1.1
      params$scale <- 300
      params$beta_age <- 0.25
    },
    "Stabilized_Progressed" = {
      params$shape <- 1.2
      params$scale <- 200
      params$beta_age <- 0.15
      params$beta_trt <- 0.1
    },
    "Stabilized_Death" = {
      params$shape <- 1.0
      params$scale <- 500
      params$beta_age <- 0.3
    },
    "Progressed_Death" = {
      params$shape <- 1.8
      params$scale <- 120
      params$beta_age <- 0.25
      params$beta_bmi <- 0.1
    }
  )

  params
}
