# tsgarch developer notes

## Verification

- Fast iteration on R-only changes: `devtools::load_all()` then
  `testthat::test_file("tests/testthat/test-<name>.R")`. The compiled TMB and
  Rcpp objects in `src/` are checked in, so no rebuild is needed unless a
  `.cpp`/`.hpp` file changed.
- Changing anything under `src/TMB/` requires recompiling
  `tsgarch_TMBExports` (`R CMD INSTALL .`, or `src/TMB/compile.R`). This is
  slow (the resulting object file is ~85MB) - avoid triggering it needlessly.
- **`.hpp` edits are silently skipped by the build.** `TMB::compile()` only
  compares the mtime of `tsgarch_TMBExports.cpp` against the object file, and
  the model templates are `.hpp` files included by it. Editing a template and
  running `R CMD INSTALL .` therefore leaves a stale DLL in place, which
  surfaces as a confusing runtime error (e.g. "Names in map must correspond
  to parameter names" when a new `PARAMETER_VECTOR` was added). Always
  `touch src/TMB/tsgarch_TMBExports.cpp` before recompiling after a `.hpp`
  change.
- Test files worth running as regression guards for any change to the
  estimation/filtering/simulation core: `test-estimation.R`,
  `test-filtering.R`, `test-simulation.R`, `test-arma.R`.
- `tests/longtests/` holds the slow backtest/prediction/profile tests; they
  are not run by `R CMD check`.

## Parameter naming conventions

Parameters live in a model's `parmatrix` (`parameter`, `value`, `lower`,
`upper`, `estimate`, `scale`, `group`, `equation`, `symbol`). Greek symbols
already in use, and what they mean, so that new parameters do not collide:

| symbol | group / parameter | meaning |
|---|---|---|
| `\mu` | `mu` | unconditional mean |
| `r^{ar}`, `r^{ma}` | `arpacf`, `mapacf` | raw Durbin-Levinson (pacf) ARMA parameters |
| `\tau` | `tau` | mean equation (`xreg`) regressor coefficients |
| `\omega` | `omega` | variance intercept |
| `\alpha`, `\beta` | `alpha`, `beta` | ARCH / GARCH |
| `\gamma`, `\eta`, `\delta` | `gamma`, `eta`, `delta` | asymmetry / rotation / power terms |
| `\rho`, `\phi` | `rho`, `phi` | cgarch permanent / transitory component |
| `\xi` | `xi` | variance equation (`vreg`) regressor coefficients |
| `\zeta`, `\nu`, `\lambda` | `distribution` | skew, shape, lambda |

Also reserved, though not parameters: `\chi_{k,t}` denotes the variance
regressor *data* throughout `vignettes/garch_models.Rmd`, `\kappa` is the
`ADREPORT`ed expected-absolute-moment term of the asymmetric models, and
`\psi` denotes the ARMA impulse response weights in `arma_irf()`.

Note for cross-package consistency: the sibling package **tsarma** names its
*mean* equation regressor coefficients `\xi`, which in tsgarch is already
taken by the *variance* equation regressors. tsgarch therefore uses `\tau`
for mean regressors; tsarma's `xi` is tsgarch's `tau`.

Ordering matters: `pscale` is handed to the TMB templates as
`parmatrix$scale` in row order, and each template walks it positionally by
group. A new parameter group must be inserted in the parmatrix at the same
position at which the templates read it.

## Regressor validation policy (mean equation `xreg` / `newxreg`)

At specification time (`check_xreg()`): number of rows must match `y`, the
time index must match, and no `NA`/`NaN`/`Inf` values are accepted. A
rank-deficient `xreg` (including a constant column when `constant = TRUE`,
which is collinear with `mu`) is rejected.

At prediction/simulation/filtering time:

- `NCOL(newxreg)` must equal `NCOL(xreg)` as supplied at estimation - error
  otherwise.
- the number of rows must match the horizon exactly (`h` for
  `predict()`/`simulate()`, `NROW(y)` for `tsfilter()`) - error otherwise.
- no `NA`/`NaN`/`Inf` values - error otherwise.
- if the model was estimated with `xreg` but `newxreg` is `NULL`, a zero
  matrix is substituted and a warning is emitted (rather than erroring). This
  matches the existing `newvreg` behaviour of `tsfilter()`
  (`.check_y_filter()`), and sits between rugarch (silent zero-padding, see
  `.forcregressors()`) and tsarma (hard error).

Note the resulting asymmetry: the *variance* regressor path
(`.process_prediction_regressors()` with `newvreg`) still errors on a missing
argument in `predict()`. Aligning `newvreg` with the warn-and-zero policy is
a pending decision, not an oversight.

Note also that `simulate()`'s legacy `vreg` argument is a *pre-multiplied*
numeric vector of length `h` (i.e. `v %*% xi`), whereas the mean equation
`xreg` argument is a matrix/xts of `h` rows by `NCOL(xreg)` columns which is
pre-multiplied internally. The latter follows tsarma and `predict()`'s
`newxreg`, and unlike the former it can actually be validated.

## Parked discussions

- `.spec2newspec()` (`R/utilities.R`) rebuilds a spec using
  `object$model$model`, which is `"igarch"` for a model the user specified as
  `"ewma"` (the coercion happens in `garch_modelspec()`). Values are copied,
  so `simulate()`/`tsfilter()` results are numerically identical, but the
  `estimate` flag on `omega` differs, which shifts `npars` by one and hence
  `AIC`/`BIC` for a filtered ewma model. `spec$model$model_name` now retains
  the original name and would fix this in one line, but it changes existing
  output, so it is deliberately left alone pending discussion.
- `check_newxreg()` (`R/utilities.R`) is defined but never called.
