---
title: 'panelecon: Comprehensive Panel Econometric Methods'
tags:
  - R
  - econometrics
  - panel data
  - unit root tests
  - cointegration
  - quantile regression
  - cross-sectional dependence
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 25 March 2026
bibliography: paper.bib
---

# Summary

`panelecon` is a comprehensive R toolkit for panel data econometrics, consolidating a broad set of methods that are typically scattered across multiple Stata ado-files or specialized software. The package covers the full pipeline of modern panel analysis: pre-testing for unit roots and parameter heterogeneity (Hsiao homogeneity tests, Swamy parameter heterogeneity test); cross-sectional dependence testing in quantile regressions; panel Granger causality tests including Fourier Toda-Yamamoto and quantile Granger causality; group-mean and pooled fully modified OLS (FMOLS) estimation; quantile regression slope homogeneity tests; panel quantile unit root tests with common shocks and structural breaks; fixed effects vector decomposition and filtered estimators; and missing data detection and imputation tools for unbalanced panels.

# Statement of Need

Modern panel econometrics demands rigorous pre-testing for cross-sectional dependence, slope heterogeneity, and non-stationarity before any inference step. While packages such as `plm` [@Croissant2008] handle standard linear panel models, advanced tests like Fourier Toda-Yamamoto Granger causality, quantile panel unit roots with structural breaks, and quantile common correlated effects are unavailable in R. Researchers studying macroeconomic convergence, energy-growth nexuses, and financial integration across countries have had to rely on Stata commands with limited scriptability and reproducibility. `panelecon` provides a reproducible, script-friendly R alternative covering the pre-testing, estimation, and inference stages of panel analysis within a single package.

# Usage

## Pre-Testing: Heterogeneity and Cross-Sectional Dependence

```r
library(panelecon)

# Swamy parameter heterogeneity test and Hsiao homogeneity tests
result_pre <- xtpretest(y ~ x1 + x2, data = panel_df,
                        index = c("id", "time"))
print(result_pre)

# Cross-sectional dependence in quantile regression
result_csd <- xtcsdq(y ~ x1 + x2, data = panel_df,
                     index = c("id", "time"), tau = 0.5)
print(result_csd)
```

## Panel Granger Causality

```r
# Fourier Toda-Yamamoto and quantile panel Granger causality
result_caus <- xtpcaus(y ~ x, data = panel_df,
                       index = c("id", "time"), freq = 1)
summary(result_caus)
```

## Panel FMOLS Estimation

```r
# Group-mean and pooled FMOLS
result_fmols <- xtpcmg(y ~ x1 + x2, data = panel_df,
                       index = c("id", "time"), estimator = "mg")
print(result_fmols)
```

## Panel Quantile Unit Root Tests

```r
# Quantile unit root with common shocks and structural breaks
result_qur <- xtpqroot(y, data = panel_df, index = c("id", "time"),
                        tau = c(0.25, 0.5, 0.75))
summary(result_qur)
```

## Missing Data Tools

```r
# Detect and impute missing data in panel
result_mis <- xtmispanel(panel_df, index = c("id", "time"),
                          method = "linear")
```

# Implementation

`panelecon` is written in pure R with optional dependencies on `plm` and `zoo` for panel data structures. The `xtpretest()` function implements the Swamy [-@Swamy1970] statistic alongside Hsiao-type homogeneity tests. The `xtcsdq()` function extends the Pesaran [-@Pesaran2004] CD test to quantile regression residuals. Fourier Toda-Yamamoto causality testing in `xtpcaus()` augments the VAR with trigonometric terms to capture smooth structural shifts. The `xtpcmg()` function implements mean-group and pooled FMOLS following Pedroni [-@Pedroni2000]. Panel quantile unit root tests in `xtpqroot()` follow the framework accommodating common factors and structural breaks. The `xtqsh()` function implements quantile slope homogeneity tests. Fixed effects vector decomposition is implemented in line with Plümper and Troeger [-@Plumper2007].

# References
