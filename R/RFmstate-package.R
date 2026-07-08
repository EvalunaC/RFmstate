#' @title RFmstate: Random Forest-Based Multistate Survival Analysis
#'
#' @description Fits cause-specific random survival forests for flexible
#'   multistate survival analysis with covariate-adjusted transition
#'   probabilities computed via product-integral. For each transient state,
#'   competing transitions are modeled by separate random forests, and
#'   patient-specific transition probability matrices are assembled from
#'   the predicted cumulative hazards using the product-integral formula.
#'   Also provides a standalone Aalen-Johansen nonparametric estimator as
#'   a covariate-free baseline. Supports arbitrary state spaces with any
#'   number of states (three or more) and any set of allowed transitions,
#'   applicable to clinical trials, disease progression, reliability
#'   engineering, and other domains where subjects move among discrete
#'   states over time. The package provides:
#'   \itemize{
#'     \item State space and transition structure definition
#'     \item Wide-to-long data conversion for multistate counting processes
#'     \item Cause-specific random forest fitting per origin state
#'     \item Transition probability matrices via product-integral of predicted
#'       cumulative hazards
#'     \item Aalen-Johansen nonparametric estimation (covariate-free baseline)
#'     \item Per-transition feature importance
#'     \item Bias-variance diagnostics with Brier score and C-index
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
