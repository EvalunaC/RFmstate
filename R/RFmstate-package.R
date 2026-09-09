#' @title RFmstate: Random Forest-Based Multistate Survival Analysis
#'
#' @description Fits transition-specific cause-specific random survival forests
#'   on a clock-reset duration scale for acyclic, non-recurrent multistate
#'   processes. Patient/profile state probabilities are assembled by
#'   semi-Markov convolution and are conditional on fresh state entry. The
#'   package supports a common initial state, one recorded entry per state,
#'   baseline time-fixed covariates, competing exits, and independent right
#'   censoring. It does not support left truncation, cycles/recurrent visits,
#'   time-dependent covariates, or ongoing-sojourn dynamic prediction. A
#'   strict predictor contract prevents IDs, event/censoring times, response
#'   fields, and arbitrary long-format columns from entering a forest. Every
#'   full-data or cross-validation fit learns its predictor schema only from
#'   its own fitting rows, and every reported OOB statistic has verified OOB
#'   coverage. A calendar-time Aalen-Johansen point estimator is provided as a
#'   covariate-free descriptive baseline. The package provides:
#'   \itemize{
#'     \item State space and transition structure definition
#'     \item Wide-to-long data conversion for multistate counting processes
#'     \item Cause-specific random forest fitting per origin state
#'     \item Entry-conditioned state probabilities via semi-Markov convolution
#'     \item Aalen-Johansen point estimation (covariate-free baseline)
#'     \item Per-transition feature importance
#'     \item Genuine edge OOB concordance and patient-level cross-validated
#'       IPCW Brier scores
#'     \item Comprehensive visualizations
#'   }
#'
#' @docType package
#' @name RFmstate-package
#' @aliases RFmstate
#'
#' @importFrom survival Surv
#' @importFrom ranger ranger
#' @importFrom stats predict quantile var median pweibull rweibull rbinom
#'   runif rnorm complete.cases model.matrix
#' @importFrom graphics plot lines legend par axis mtext polygon abline
#'   barplot text box layout arrows rect title image points
#' @importFrom grDevices adjustcolor rgb
#' @importFrom utils head tail
"_PACKAGE"
