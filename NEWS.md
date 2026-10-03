# panelecon 1.1.0

This release synchronises every component with the audited standalone
package that holds the same functions. The standalone packages were checked
against their source papers and corrected; panelecon 1.0.1 still carried
the older code for two of them.

## Replaced components

* `xtpcmg()` now matches caustests 1.1.4 (previously the code of caustests
  1.1.3). Corrections taken over from caustests 1.1.4:
  * The one-sided long-run covariance used in the fully modified bias
    correction was transposed (it estimated the sum of E(u_t v_{t+j})
    instead of E(v_t u_{t+j})), and the quadratic spectral and Daniell
    kernels did not use the one-sided weights of the authors' code;
    group-mean and pooled FM-OLS estimates were therefore biased. The
    long-run covariance now follows the authors' lr_varmod.m.
  * Pooled model: the covariance matrix is now the asymptotic covariance of
    de Jong and Wagner (2022) for one-way and two-way effects (as in the
    authors' PanelEKC code), with a heteroskedasticity-robust sandwich for
    the controls; the previous version used sigma^2 (X'X)^-1 from the FM
    residuals, and a unit matrix when X'X was singular.
  * Cross-section robust covariance (`corr_rob = TRUE`): uses the
    conditional long-run covariance between units instead of the
    covariance of u alone.
  * All `xtpcmg()` estimates and standard errors now reproduce the Stata
    module xtpcmg 1.0.2; the corresponding test was added
    (`tests/testthat/test-xtpcmg-stata.R`). The function signature is
    unchanged.
* `xtfifevd()` now matches xtfifevd 1.1.0. panelecon 1.0.1 contained an
  independent single-file implementation of FEVD, FEF and FEF-IV that did
  not correspond to any release of the standalone package; it has been
  replaced by the audited code (files `R/xtfifevd_xtfifevd.R`,
  `R/xtfifevd_estimation.R`, `R/xtfifevd_methods.R`,
  `R/xtfifevd_diagnostics.R`). Corrections taken over from xtfifevd 1.1.0:
  * Formula handling: the right-hand side is split at the top-level `|` and
    each part is built with `model.frame()` and `model.matrix()`, so
    transformations such as `log(x)`, `I(x^2)`, `x:w` and `log(y)` and
    factors (with their contrasts) are handled as in `lm()`. Duplicated
    `(id, time)` pairs are an error.
  * Variance of the time-varying coefficients `beta`: the panel-robust
    matrix of Pesaran and Zhou (2018, eq. 18) is reported by default
    (`vcov_beta = "robust"`); `vcov_beta = "classical"` gives
    `sigma2_e (X'MX)^-1`. The chosen matrix is also used inside the Pesaran
    and Zhou variance of `gamma` and of the intercept.
  * Full covariance matrix: `vcov()` contains the covariance between
    `gamma` and `beta` from Pesaran and Zhou (2018, eq. A.11) and the
    variance of and covariances with the intercept by the delta method;
    the same derivation is applied to FEF-IV with the instrument
    projections (eq. 48 and 51). The previous intercept standard errors
    (stage 2 OLS for FEF, 2SLS for FEF-IV) are no longer used.
  * `fevd()` runs the three stages of Plumper and Troeger (2007) literally;
    the stage 3 coefficient on the unexplained unit effect is returned in
    `delta` (equal to 1 by construction) and the naive stage 3 standard
    errors in `stage3$se_naive` for reference only. Inference uses the
    Pesaran and Zhou standard errors.
  * Rarely changing variables after `|`: the unit (panel) mean is used in
    all stages, with a warning.
  * `sigma2_u` is the variance of the unexplained unit effect (stage 2
    residual variance minus `sigma2_e * mean(1 / T_i)`, truncated at zero).
  * New exported functions `fevd()`, `fef()`, `fef_iv()` and `bw_ratio()`,
    and new S3 methods `coef()`, `vcov()`, `confint()` and
    `print.summary.xtfifevd()`. `bw_ratio()` reports between/within SD
    ratios and prints the thresholds of Plumper and Troeger (2007, Fig. 4)
    for orientation only.
  * Breaking change in the interface of `xtfifevd()`. The old call
    `xtfifevd(y ~ x, data, index = c("id", "time"), zinvariants = "z",
    method, instruments = c("r1"), robust)` becomes
    `xtfifevd(y ~ x | z, data, id = "id", time = "time", method,
    instruments = ~ r1, vcov_beta)`. `method = "fefiv"` is now
    `"fef_iv"`. The returned object has elements `coefficients`, `vcov`,
    `beta`, `gamma`, `intercept`, `V_gamma_pz`, `sigma2_e`, `sigma2_u`,
    `N`, `N_g`, `T_bar` and others (see `?xtfifevd`) instead of `beta_fe`,
    `se_beta`, `se_gamma`, `se_alpha`, `alpha`, `T_avg`. `summary()` now
    returns a `"summary.xtfifevd"` object with a coefficient table.
  * The tests of xtfifevd 1.1.0 were added
    (`tests/testthat/test-xtfifevd.R`,
    `tests/testthat/test-xtfifevd-pesaran-zhou.R`); they pin the FEF and
    FEF-IV coefficients and the full covariance matrix to hand computations
    of the Pesaran and Zhou formulas and `delta = 1`.

## Components verified identical to their standalone source

* `xtcsdq()`, `xtmispanel()`, `xtpretest()` and `xtqsh()`: identical to
  paneltests 1.0.6.
* `xtpcaus()`: identical to caustests 1.1.4 (includes the Toda-Yamamoto lag
  selection and Fourier index corrections of caustests 1.1.3).
* `xtpqroot()`: identical to xtpqroot 1.0.0 (code); only the example
  wrapping differs.

## Other changes

* Documentation: all Rd files are now generated by roxygen2 from the
  comments in `R/`. The `\value` section of `?xtpqroot` was missing in the
  generated file because of unescaped percent signs; the three bootstrap
  critical value entries now say "percent". The `\describe` blocks of
  `?xtpcaus` and `?xtpqroot` no longer contain text between items.
* "and" replaces "&" between author names in the documentation.
* New imports from stats (`as.formula`, `confint`, `fitted`,
  `printCoefmat`, `residuals`, `terms`) and utils (`packageVersion`) for
  the xtfifevd component.
