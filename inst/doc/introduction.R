## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

# One configuration is used throughout the vignette. These values keep the
# source build practical; they are not inferential adequacy recommendations.
canonical_seed <- 42L
canonical_n <- 300L
canonical_trees <- 200L
canonical_min_events <- 3L
canonical_sparse_warning <- 20L
canonical_covariates <- c("age", "sex", "BMI", "treatment")

## ----define-states------------------------------------------------------------
library(RFmstate)

# Use the built-in clinical trial structure
ms <- clinical_states()
print(ms)

## ----custom-states, eval=FALSE------------------------------------------------
# # A simple 3-state illness-death model
# ms_simple <- define_multistate(
#   state_names = c("Healthy", "Sick", "Dead"),
#   absorbing = "Dead",
#   transitions = list(
#     Healthy = c("Sick", "Dead"),
#     Sick = c("Dead")
#   )
# )
# 
# # A 4-state model with recovery
# ms_recovery <- define_multistate(
#   state_names = c("Healthy", "Sick", "Recovered", "Dead"),
#   absorbing = "Dead",
#   transitions = list(
#     Healthy = c("Sick", "Dead"),
#     Sick = c("Recovered", "Dead"),
#     Recovered = c("Dead")
#   )
# )

## ----simulate-----------------------------------------------------------------
dat <- sim_clinical_data(
  n = canonical_n, structure = ms, seed = canonical_seed
)
head(dat)

## ----prepare------------------------------------------------------------------
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
  covariates = canonical_covariates
)
print(msdata)
head(msdata)

## ----aj-----------------------------------------------------------------------
aj <- aalen_johansen(msdata)
print(aj)

## ----aj-plot, fig.cap="State occupation probabilities from Aalen-Johansen estimator"----
plot(aj, type = "state_occupation")

## ----aj-hazard, fig.cap="Nelson-Aalen cumulative hazards"---------------------
plot(aj, type = "cumulative_hazard")

## ----fit----------------------------------------------------------------------
fit <- rfmstate(
  msdata,
  num.trees = canonical_trees,
  min_events = canonical_min_events,
  sparse_warning = canonical_sparse_warning,
  seed = canonical_seed
)
print(fit)

## ----summary------------------------------------------------------------------
summary(fit)

## ----importance, fig.cap="Feature importance per transition"------------------
imp <- importance(fit)
print(imp)
plot(imp, type = "barplot")

## ----importance-heat, fig.cap="Feature importance heatmap"--------------------
plot(imp, type = "heatmap")

## ----predict, fig.cap="Predicted state occupation for two patient profiles"----
newdata <- data.frame(
  age = c(50, 70),
  sex = c(0, 1),
  BMI = c(24, 32),
  treatment = c(1, 0)
)

prediction_horizon <- min(fit$max_duration_by_origin)
pred <- predict(fit, newdata = newdata,
                times = seq(0, prediction_horizon, length.out = 37))

# Plot for patient 1 (young, treated)
plot(pred, type = "state_occupation", subject = 1)

# Plot for patient 2 (older, untreated)
plot(pred, type = "state_occupation", subject = 2)

## ----diagnostics--------------------------------------------------------------
diag <- diagnose(fit)
print(diag)

## ----diag-concordance, fig.cap="Concordance index per transition"-------------
plot(diag, type = "concordance")

## ----diag-cv, eval=FALSE------------------------------------------------------
# cv_diag <- diagnose(fit, method = "cv", folds = 5,
#                     eval_times = seq(0, prediction_horizon * 0.8,
#                                      length.out = 9))
# plot(cv_diag, type = "brier")

## ----diagram, fig.cap="Transition diagram with event counts"------------------
plot_transition_diagram(ms, msdata)

## ----direct-probability-------------------------------------------------------
simple_ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
elapsed_grid <- seq(0, 2, length.out = 2001)
simple_hazards <- list(
  "A->B" = data.frame(time = elapsed_grid,
                       hazard = 0.4 * elapsed_grid)
)
simple_prob <- compute_trans_prob(
  simple_hazards, simple_ms, times = c(0, 1, 2),
  target_grid_points = 512
)
simple_prob$state_occ

## ----session-info-------------------------------------------------------------
sessionInfo()

