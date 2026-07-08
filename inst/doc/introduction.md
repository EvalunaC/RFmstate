---
title: "Introduction to RFmstate"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to RFmstate}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Overview

**RFmstate** fits cause-specific random survival forests for multistate
survival analysis with covariate-adjusted transition probabilities computed
via product-integral. For each transient state, competing transitions are
modeled by separate random forests, and patient-specific transition
probability matrices are assembled from the predicted cumulative hazards
using the product-integral formula. The package also provides a standalone
Aalen-Johansen nonparametric estimator as a covariate-free baseline. This
approach is particularly suited for clinical trial data where patients
transition through discrete health states (e.g., response, progression,
death) with right censoring and competing transitions.

The package provides:

- State space and transition structure definition
- Wide-to-long data conversion for counting-process format
- Cause-specific random forest fitting per origin state
- Transition probability matrices via product-integral of predicted
  cumulative hazards
- Aalen-Johansen nonparametric estimation (covariate-free baseline)
- Per-transition feature importance
- Bias-variance diagnostics with Brier score and concordance index
- Comprehensive visualizations

## Quick Start

### 1. Define the Multistate Structure


``` r
library(RFmstate)

# Use the built-in clinical trial structure
ms <- clinical_states()
print(ms)
#> Multistate Structure
#>   States: Baseline -> Responded -> Unresponded -> Stabilized -> Progressed -> Death 
#>   Absorbing: Death 
#>   Transitions: 12 
#>     1: Baseline -> Responded
#>     2: Baseline -> Unresponded
#>     3: Baseline -> Death
#>     4: Responded -> Stabilized
#>     5: Responded -> Progressed
#>     6: Responded -> Death
#>     7: Unresponded -> Stabilized
#>     8: Unresponded -> Progressed
#>     9: Unresponded -> Death
#>     10: Stabilized -> Progressed
#>     11: Stabilized -> Death
#>     12: Progressed -> Death
```

Or define a custom structure with any number of states (3 or more):


``` r
# A simple 3-state illness-death model
ms_simple <- define_multistate(
  state_names = c("Healthy", "Sick", "Dead"),
  absorbing = "Dead",
  transitions = list(
    Healthy = c("Sick", "Dead"),
    Sick = c("Dead")
  )
)

# A 4-state model with recovery
ms_recovery <- define_multistate(
  state_names = c("Healthy", "Sick", "Recovered", "Dead"),
  absorbing = "Dead",
  transitions = list(
    Healthy = c("Sick", "Dead"),
    Sick = c("Recovered", "Dead"),
    Recovered = c("Dead")
  )
)
```

This pipeline works with any defined multistate structure.

### 2. Simulate Data


``` r
set.seed(42)
dat <- sim_clinical_data(n = 300, structure = ms)
head(dat)
#>   ID  age sex  BMI treatment time_Responded time_Unresponded time_Stabilized
#> 1  1 76.5   0 19.5         0           23.8               NA            27.9
#> 2  2 53.2   0 24.7         1           17.7               NA            52.7
#> 3  3 64.4   1 26.9         0             NA             16.0            46.6
#> 4  4 67.6   0 24.0         1             NA             25.4              NA
#> 5  5 64.9   1 26.5         0             NA             42.2              NA
#> 6  6 58.7   0 24.4         1            6.1               NA             9.0
#>   time_Progressed time_Death time_censored
#> 1           274.8      275.7            NA
#> 2           164.3      225.9            NA
#> 3            73.7      113.6            NA
#> 4            62.6      154.1            NA
#> 5              NA       61.2            NA
#> 6           108.4      121.0            NA
```

### 3. Prepare Multistate Data

Convert wide-format data to long format:


``` r
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
print(msdata)
#> Multistate Data (msdata)
#>   Patients: 300 
#>   Intervals: 936 
#>   Transitions observed: 910 
#>   Censored intervals: 26 
#>   States: Baseline, Responded, Unresponded, Stabilized, Progressed, Death 
#> 
#> Transition counts:
#>              to
#> from          Death Progressed Responded Stabilized Unresponded
#>   Baseline       23          0       193          0          84
#>   Progressed    172          0         0          0           0
#>   Responded      17         36         0        140           0
#>   Stabilized     50        111         0          0           0
#>   Unresponded    12         31         0         41           0
```

### 4. Aalen-Johansen Nonparametric Baseline

Compute nonparametric (covariate-free) transition probability estimates:


``` r
aj <- aalen_johansen(msdata)
print(aj)
#> Aalen-Johansen Estimate
#>   Time range: [0.3, 358.2]
#>   Event times: 728 
#>   States: Baseline, Responded, Unresponded, Stabilized, Progressed, Death 
#> 
#> Event counts per transition:
#>         from          to n_events
#>     Baseline   Responded      193
#>     Baseline Unresponded       84
#>     Baseline       Death       23
#>    Responded  Stabilized      140
#>    Responded  Progressed       36
#>    Responded       Death       17
#>  Unresponded  Stabilized       41
#>  Unresponded  Progressed       31
#>  Unresponded       Death       12
#>   Stabilized  Progressed      111
#>   Stabilized       Death       50
#>   Progressed       Death      172
#> 
#> Final state occupation probabilities:
#>   Baseline: 0
#>   Responded: 0
#>   Unresponded: 0
#>   Stabilized: 0
#>   Progressed: 0.0138
#>   Death: 0.9862
```


``` r
plot(aj, type = "state_occupation")
```

![State occupation probabilities from Aalen-Johansen estimator](figure/aj-plot-1.png)


``` r
plot(aj, type = "cumulative_hazard")
```

![Nelson-Aalen cumulative hazards](figure/aj-hazard-1.png)

### 5. Fit Random Forest Model


``` r
fit <- rfmstate(
  msdata,
  covariates = c("age", "sex", "BMI", "treatment"),
  num.trees = 200,
  seed = 42
)
print(fit)
#> Random Forest Multistate Model
#> Call: rfmstate(msdata = msdata, covariates = c("age", "sex", "BMI", 
#>     "treatment"), num.trees = 200, seed = 42)
#> 
#> Covariates: age, sex, BMI, treatment 
#> Parameters:
#>   num.trees: 200 
#>   mtry: 2 
#>   min.node.size: 15 
#> 
#> Models fitted per origin state:
#>   Baseline (n=300): -> Responded, Unresponded, Death
#>   Responded (n=193): -> Stabilized, Progressed, Death
#>   Unresponded (n=83): -> Stabilized, Progressed, Death
#>   Stabilized (n=181): -> Progressed, Death
#>   Progressed (n=178): -> Death
```

### 6. Model Summary


``` r
s <- summary(fit)
```

### 7. Feature Importance


``` r
imp <- importance(fit)
print(imp)
#> Feature Importance per Transition
#> ============================================================ 
#> 
#>           Baseline -> Responded Baseline -> Unresponded Baseline -> Death
#> age                      0.0155                 -0.0039            0.0574
#> sex                      0.0013                 -0.0024           -0.0103
#> BMI                     -0.0003                  0.0318           -0.0077
#> treatment                0.0354                  0.0356           -0.0088
#>           Responded -> Stabilized Responded -> Progressed Responded -> Death
#> age                        0.0019                 -0.0182             0.0002
#> sex                       -0.0004                 -0.0129            -0.0020
#> BMI                        0.0081                 -0.0009             0.0300
#> treatment                  0.0064                 -0.0091             0.0524
#>           Unresponded -> Stabilized Unresponded -> Progressed
#> age                         -0.0014                    0.0134
#> sex                         -0.0118                   -0.0106
#> BMI                         -0.0024                    0.0200
#> treatment                   -0.0134                    0.0008
#>           Unresponded -> Death Stabilized -> Progressed Stabilized -> Death
#> age                        NaN                   0.0440              0.0349
#> sex                        NaN                   0.0054              0.0141
#> BMI                        NaN                  -0.0045             -0.0088
#> treatment                  NaN                   0.0095              0.0155
#>           Progressed -> Death
#> age                    0.0778
#> sex                   -0.0051
#> BMI                   -0.0036
#> treatment             -0.0016
#> 
#> Top variables per transition:
#>   Baseline -> Responded: treatment (0.0354)
#>   Baseline -> Unresponded: treatment (0.0356)
#>   Baseline -> Death: age (0.0574)
#>   Responded -> Stabilized: BMI (0.0081)
#>   Responded -> Progressed: BMI (-9e-04)
#>   Responded -> Death: treatment (0.0524)
#>   Unresponded -> Stabilized: age (-0.0014)
#>   Unresponded -> Progressed: BMI (0.02)
#>   Stabilized -> Progressed: age (0.044)
#>   Stabilized -> Death: age (0.0349)
#>   Progressed -> Death: age (0.0778)
plot(imp, type = "barplot")
```

![Feature importance per transition](figure/importance-1.png)


``` r
plot(imp, type = "heatmap")
```

![Feature importance heatmap](figure/importance-heat-1.png)

### 8. Predict for New Patients


``` r
newdata <- data.frame(
  age = c(50, 70),
  sex = c(0, 1),
  BMI = c(24, 32),
  treatment = c(1, 0)
)

pred <- predict(fit, newdata = newdata, times = seq(10, 365, by = 10))

# Plot for patient 1 (young, treated)
plot(pred, type = "state_occupation", subject = 1)
```

![Predicted state occupation for two patient profiles](figure/predict-1.png)

``` r

# Plot for patient 2 (older, untreated)
plot(pred, type = "state_occupation", subject = 2)
```

![Predicted state occupation for two patient profiles](figure/predict-2.png)

### 9. Diagnostics


``` r
diag <- diagnose(fit)
print(diag)
#> RF Multistate Model Diagnostics
#> ============================================================ 
#> 
#> OOB Prediction Error:
#> ---------------------------------------- 
#>   Baseline -> Responded     0.4220
#>   Baseline -> Unresponded   0.3877
#>   Baseline -> Death         0.4372
#>   Responded -> Stabilized   0.4698
#>   Responded -> Progressed   0.5251
#>   Responded -> Death        0.3142
#>   Unresponded -> Stabilized 0.5913
#>   Unresponded -> Progressed 0.5635
#>   Unresponded -> Death      0.6269
#>   Stabilized -> Progressed  0.4190
#>   Stabilized -> Death       0.4940
#>   Progressed -> Death       0.3837
#> 
#> Concordance Index (C-index):
#> ---------------------------------------- 
#>   Baseline -> Responded     0.7660
#>   Baseline -> Unresponded   0.8690
#>   Baseline -> Death         0.9219
#>   Responded -> Stabilized   0.7751
#>   Responded -> Progressed   0.8911
#>   Responded -> Death        0.9333
#>   Unresponded -> Stabilized 0.8242
#>   Unresponded -> Progressed 0.8612
#>   Unresponded -> Death      0.8916
#>   Stabilized -> Progressed  0.8187
#>   Stabilized -> Death       0.8879
#>   Progressed -> Death       0.7800
#> 
#> Bias-Variance Decomposition:
#> ------------------------------------------------------------ 
#>   Transition                    Bias      Var      MSE
#> ------------------------------------------------------------ 
#>   Baseline -> Responded       0.0204   0.0295   0.1634
#>   Baseline -> Unresponded     0.0088   0.0185   0.1022
#>   Baseline -> Death           0.0149   0.0084   0.0374
#>   Responded -> Stabilized     0.0009   0.0242   0.1657
#>   Responded -> Progressed     0.0340   0.0104   0.0896
#>   Responded -> Death          0.0075   0.0094   0.0385
#>   Unresponded -> Stabilized   0.0257   0.0152   0.1501
#>   Unresponded -> Progressed   0.0520   0.0132   0.1355
#>   Unresponded -> Death        0.0331   0.0058   0.0750
#>   Stabilized -> Progressed    0.0214   0.0265   0.1485
#>   Stabilized -> Death         0.0213   0.0141   0.0979
#>   Progressed -> Death        -0.0073   0.0350   0.1615
```


``` r
plot(diag, type = "brier")
```

![Time-dependent Brier score](figure/diag-brier-1.png)


``` r
plot(diag, type = "concordance")
```

![Concordance index per transition](figure/diag-concordance-1.png)


``` r
plot(diag, type = "bias_variance")
```

![Bias-variance decomposition](figure/diag-bv-1.png)

### 10. Transition Diagram


``` r
plot_transition_diagram(ms, msdata)
```

![Transition diagram with event counts](figure/diagram-1.png)

## Methodology

### Product-Integral Framework

Both the nonparametric (Aalen-Johansen) and the random forest methods in
this package compute transition probability matrices $P(s,t)$ using the
product-integral:

$$P(s,t) = \prod_{s < u \leq t} (I + dA(u))$$

where $dA(u)$ is a matrix of hazard increments at time $u$, with off-diagonal
entries $dA_{hj}(u)$ representing the $h \to j$ transition hazard increment
and diagonal entries $dA_{hh}(u) = -\sum_{j \neq h} dA_{hj}(u)$. The two
methods differ in how the hazard increments $dA(u)$ are estimated.

### Aalen-Johansen Estimator (Nonparametric Baseline)

The Aalen-Johansen (AJ) estimator is the nonparametric generalization of the
Kaplan-Meier estimator to multistate models under the Markov assumption. It
estimates hazard increments from the data directly via the Nelson-Aalen
formula:

$$d\hat{A}_{hj}(u) = \frac{dN_{hj}(u)}{Y_h(u)}$$

where $dN_{hj}(u)$ counts the observed $h \to j$ transitions at time $u$ and
$Y_h(u)$ is the number at risk in state $h$ just before time $u$. This
provides population-level transition probabilities without covariate
adjustment and serves as a covariate-free baseline in the package.

### Random Forest Multistate Approach

For covariate-adjusted predictions, we decompose the multistate model into
per-origin-state competing risks problems:

1. **For each transient state** $h$, identify all outgoing transitions
2. **Fit a cause-specific RSF**: For each destination state $j$, fit a
   random survival forest treating transition $h \to j$ as the event of
   interest and all other transitions as censored
3. **Extract cumulative hazards**: Convert each forest's predicted survival
   function $\hat{S}_{hj}(t|\mathbf{x})$ to a cumulative hazard via
   $\hat{H}_{hj}(t|\mathbf{x}) = -\log \hat{S}_{hj}(t|\mathbf{x})$
4. **Assemble via product-integral**: Compute hazard increments from the
   predicted cumulative hazards and apply the product-integral formula to
   obtain the full transition probability matrix $P(s,t|\mathbf{x})$

This approach leverages the flexibility of random forests to capture
nonlinear covariate effects and interactions while maintaining the
interpretability of the multistate framework. The product-integral step
ensures that the resulting transition probability matrices are coherent
(rows sum to one, non-negative entries).

### Diagnostics

- **OOB Error**: Out-of-bag prediction error from the random forest ensemble
- **C-index**: Concordance between predicted risk and observed event ordering
- **Brier Score**: Calibrated prediction error combining bias and variance
- **Bias-Variance Decomposition**: Systematic vs. random components of
  prediction error

## References

- Aalen, O.O. & Johansen, S. (1978). An empirical transition matrix for
  non-homogeneous Markov chains based on censored observations.
  *Scandinavian Journal of Statistics*, 5(3), 141-150.
- Ishwaran, H. et al. (2008). Random survival forests.
  *Annals of Applied Statistics*, 2(3), 841-860.
- Putter, H., Fiocco, M. & Geskus, R.B. (2007). Tutorial in biostatistics:
  Competing risks and multi-state models. *Statistics in Medicine*, 26,
  2389-2430.
