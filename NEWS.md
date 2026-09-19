# tsgarch 1.0.5

* Added a jointly estimated ARMA(p,q)-X mean equation, via new `arma`, `xreg`
and `xreg_type` arguments to `garch_modelspec` (defaulting to `c(0,0)` and
`NULL` respectively, so fully backward compatible). Available for all 8
GARCH flavors, including `igarch` and `ewma`. The `constant` argument
remains independent of both the ARMA order and the regressors, and continues
to control only whether the unconditional mean `mu` is estimated or fixed at
zero.
* Regressors enter the conditional mean through `xreg` (an xts matrix aligned
to `y`), with `xreg_type` selecting between two conventions: `xreg_type =
"arma_errors"` (the default) runs the ARMA recursion on
`y_t - mu - x_t'tau`, so `tau` is the long-run marginal effect - this
matches `stats::arima`'s `xreg` semantics and the sibling **tsarma**
package; `xreg_type = "armax"` (the **rugarch** convention) adds
`x_t'tau` to the conditional mean at time `t` only, so `tau` is the
impact effect and the long-run effect is `tau/(1 - sum(phi))`. The two
are algebraically identical whenever the AR order is zero. Both are
restrictions of the same distributed lag model and neither nests the other;
see the "ARMA Mean Equation" section of `vignettes/garch_models.Rmd` for
guidance on choosing between them. The coefficients are named `tau1`,
`tau2`, ... in the `parmatrix` (`xi` is already the variance-regressor
coefficient in this package; note that **tsarma** calls its mean-equation
coefficients `xi`, which is `tsgarch`'s `tau`). An `xreg` containing
`NA`/`NaN`/`Inf` values, or which is rank deficient (including a constant
column when `constant = TRUE`), is rejected at specification time.
* Mean regressors are supported throughout the model lifecycle:
`estimate`, `tsfilter` (via `newxreg`), `predict` (via `newxreg`, and
both parametric and bootstrap predictive distributions), `simulate`
(via `xreg`, a raw matrix with `h + burn` rows - unlike the legacy
`vreg` argument which is a pre-multiplied vector) and `garch_backtest`
(sliced per rolling window). When the model includes `xreg` but the
future regressors are not supplied, a zero matrix is substituted with a
warning; the variance-regressor argument `newvreg` in `predict()`
deliberately retains its existing hard-error behaviour. Note that future
regressors are treated as a single deterministic path shared across all
simulated paths, so regressor uncertainty is not propagated into the
predictive distribution.
* The AR/MA coefficients are guaranteed stationary/invertible by construction,
via a Durbin-Levinson/Jones partial-autocorrelation (PACF) reparameterization
of the raw optimization parameters, so no additional nonlinear constraints
are required; exact (autodiff, not finite-difference) Jacobians of the ARMA
mean equation are available throughout via TMB. The transformed AR/MA
coefficients can be extracted with the new `arma_coefficients()` function.
* `fitted()` now returns the genuinely time-varying conditional mean whenever
the mean equation has ARMA dynamics or regressors (previously always a
constant `mu`, now completing the documented "vector the size of y"
contract), and `residuals()` is now always exactly `y - fitted(object)` for
every model, including a fix for models with no ARMA (a small, purely
internal simplification with no behavior change there). A pure regression
with GARCH errors (`arma = c(0,0)` with `xreg`) therefore reports the
regression fit rather than a constant.
* `predict()`, `simulate()` and `tsfilter()` all account for the full
ARMA-X mean equation: point forecasts, simulated paths and incremental
filtering correctly reflect both the AR/MA dynamics and the regressor
contribution. The combined pre-sample/burn-in length used internally is
`max(garch order, arma order)`.
* `summary()` (and its console `print()`/`as_flextable()` methods) now
displays the transformed `ar`/`ma` coefficients with their delta-method
standard errors, rather than the raw (not directly interpretable)
Durbin-Levinson parameters used internally during estimation.
* An entire AR and/or MA polynomial in the ARMA mean equation can now be
fixed at target coefficients rather than estimated: set `value` to the
desired ar/ma coefficients (not the internal pacf parameterization) and
`estimate = 0` on every `arpacf#`/`mapacf#` row of the relevant group in the
spec's `parmatrix`; the correct internal Durbin-Levinson transform is then
applied automatically wherever needed (estimation, `arma_coefficients()`,
`predict()`, `simulate()`, `tsfilter()`), and the fixed coefficients are
validated for stationarity/invertibility up front. Fixing only some (not
all) of the lags of a given polynomial is not supported, since the
reparameterization couples all of a polynomial's lags together, and raises
an informative error, as does mixing `estimate` values within one group.
* Added a new ARMA diagnostic panel to `plot.tsgarch.estimate()`
(`plot(object, type = "arma")`), for models with a non-zero `arma` order.
The panel shows: inverse AR/MA roots against the unit circle (with visual
flags for near common-factor cancellation), the impulse response function
(plotted from lag 1, since `psi_0 = 1` identically and carries no
information about the fitted model; the optional `cumulative` line is
accumulated from lag 1 as well, so both stay on the same scale),
and the ACFs of the standardized residuals `z_t` and `z_t^2`. All four
panels are shown by default (`which = NULL`); individual panels can be
selected via `which`. The residual ACF panels support three envelopes:
the default asymptotic `"bartlett"` band, a `"simulate"` band from the
fitted innovation distribution at the point parameter estimates (no
parameter uncertainty), and a `"parametric"` band that perturbs the full
parameter vector using `vcov()` and re-filters the observed data to reflect
parameter estimation uncertainty. `type = "garch"` (the original
volatility/news-impact/QQ panel) remains the default and is unchanged. New
computational helpers `arma_inverse_roots()`, `arma_irf()` and
`arma_near_cancellation()` are exported for programmatic use.
* Fixed `tsbacktest()` silently refitting the wrong model inside its
rolling windows: the inner `garch_modelspec()` call previously passed
neither `model` nor `arma`, so e.g. an `egarch` spec was backtested as a
vanilla `garch` and any ARMA mean equation was dropped. The spec now
records the user-facing model name in `model$model_name` (necessary
because `ewma` is coerced to `igarch` internally) and the backtest refits
each window with the original model and ARMA order.
* Fixed `predict(..., sim_method = "bootstrap")` ignoring the ARMA mean
equation: `garch_bootstrap()` now uses the combined pre-sample length
`max(garch order, arma order)` and passes `series_init`/`resid_init` to
`simulate()`, so the bootstrap predictive distribution is centered on the
ARMA mean forecast rather than on the unconditional mean `mu`.
* Fixed `constant_variance` under `variance_targeting = TRUE` being
computed as the mean squared deviation from the constant `mu` instead of
the unconditional variance of the ARMA innovations `eps = y -
conditional_mu` (see the "Variance Targeting" section of
`vignettes/garch_models.Rmd`); the R-side `unconditional()`/`target_omega`
now agree with the TMB-side variance target when `arma != c(0,0)`. Both
the `estimate()` and `tsfilter()` code paths were fixed; non-ARMA models
are numerically unchanged.
* Fixed `simulate()` for `igarch`/`ewma` models with `arma != c(0,0)`
dropping the ARMA mean equation entirely - the simulated series was
just `mu + eps`. The ARMA overlay is now applied, which also corrects
the simulated predictive distribution returned by
`predict(..., nsim > 0)` for these models.
* Fixed `tsfilter()`'s zero fill for a missing `newvreg` being built
with `as.matrix(0, ...)` (a 1x1 matrix) instead of
`matrix(0, nrow, ncol)`.
* Fixed `simulate()` on a spec taken from an estimated object
(`mod$spec`, which carries `parmatrix = NULL`) aborting the R session;
it now raises a clean error asking the user to assign the estimated
parmatrix onto the spec first.
* `check_xreg()` now rejects `Inf` values (previously only `NA`/`NaN`),
and its time-index match check, which previously never fired because
`all.equal()` returns a description string rather than `FALSE`, now
works as intended.
* Fixed spec re-specification dropping the `ewma` restriction:
`garch_modelspec` coerces `model = "ewma"` to `igarch` internally (the
difference is only the fixed `omega` parmatrix row), and both
`.spec2newspec()` and `garch_profile()` rebuilt specs from the coerced
name, so a filtered `ewma` model carried a spurious free `omega` (npars
and AIC/BIC off by one parameter's worth) and `garch_profile` re-fit an
igarch model. Both now re-specify from the retained user-facing
`model$model_name`. Filtered `sigma` is unchanged; filtered `ewma`
AIC/BIC shift accordingly.
* Fixed `tsprofile()` profiling the wrong mean equation: its inner
`garch_modelspec()` call passed no `arma`, so the data were simulated from
the full ARMA model while every re-estimation fitted a constant-mean
model, and the ARMA parameters never appeared in the profile at all.
Regressors in the mean equation are now refused explicitly (as `vreg`
already was), since the simulation step does not carry them into the
simulated sample.
* `simulate()` now raises an error instead of silently returning zeros when
the implied initial variance is not positive and finite. For an `ewma`
specification `omega` is fixed at zero, so the seed
`omega/(1 - 0.999)` was zero and the variance recursion stayed at zero for
every step, returning `sigma` identically zero and a series equal to `mu`.
Supply `var_init` for such models; `predict()` always did, and every other
flavour has a positive implied seed and is unaffected.
* `tsequation()` (and the `as_flextable` summary footer built from it)
now renders the conditional mean equation `eq_mean` first in the
equation block, covering the constant, ARMA terms and mean regressors
under both `xreg_type` conventions.
* A specification serialized by an older version of the package is now
handled explicitly rather than failing obscurely. `model_options` gained
the `ar`, `ma` and `xreg_type` flags during this cycle, while every
released version up to 1.0.4 wrote only its first six elements; the
padding applied on load added a single element, which left the last two
flags to be read past the end of the vector. It is padded to the full
length now. Such a specification also predates the `arpacf`, `mapacf` and
`tau` rows of `parmatrix`, and is refused with a message naming
`garch_modelspec()` as the remedy rather than the TMB error
`Error when reading the variable: 'arpacf'`.
* The hessian behind the standard errors in the scaled estimation step is
now evaluated at the solution the optimizer returns, rather than at TMB's
default of whichever point it last happened to evaluate. The two coincide
for the solver in use, so reported standard errors are unchanged.
* Fixed the pre-sample initialization of the ARCH recursion for the
asymmetric flavors (`egarch`, `gjrgarch`, `aparch`, `fgarch`). The test
deciding whether a lag lookback falls inside the pre-sample compared it
against the ARCH order rather than against the pre-sample length, and the
two coincide only when the ARCH order is the largest of the GARCH and ARMA
orders. Observations just past the pre-sample therefore took a zeroed
residual in place of the initial ARCH value, in the likelihood and in
`simulate()` alike. Both now test the pre-sample length, in the TMB
templates and in the Rcpp simulation recursions. This was already
reachable through a GARCH order with `q > p`, and the new ARMA order made
it reachable at the default `order = c(1,1)`, since the pre-sample spans
the mean recursion as well as the variance one. `logLik` and the
coefficients consequently shift for these four flavors whenever
`max(order, arma)` exceeds the ARCH order; `order = c(1,1)` with no ARMA
order is arithmetically unchanged, and `garch`, `igarch` and `ewma` are
unaffected at every order, their pre-sample ARCH input being constant. A
deterministic replication of a fitted model through `simulate()` now
recovers its `sigma` to machine precision for `gjrgarch`, `aparch` and
`fgarch` in these configurations, where the discrepancy previously reached
5e-2. `egarch` improves but remains inexact, and is still under
investigation.
* `arch_initial` is indexed by ARCH lag when passed back into `simulate()`,
so a vector shorter than the pre-sample is padded now rather than
collapsed onto its first element, which had put the lag 1 initialization
into every lag slot.
* `simulate()` now reproduces a fitted `egarch` model exactly, as the other
flavors already did. Its pre-sample branch added `alpha_j z_{t-j}` built
from the synthetic pre-sample innovations, where the likelihood contributes
nothing at all: pre-sample residuals are zeroed in the template, so the
standardized residual there is identically zero and only the `gamma_j` term
survives. The simulation matches now, which takes the discrepancy against a
fitted model's `sigma` from around 9e-3 to machine precision. Free running
`egarch` simulations shift accordingly over their first
`max(order, arma)` steps. The test that was supposed to cover this had
copied the `garch` fixtures throughout, so `egarch` simulation was in
practice never validated; it now fits an `egarch` model.
* An `egarch` simulation took its ARCH initialization from the first
`order[1]` pre-sample columns instead of all `max(order, arma)` of them.
With more than one simulated path this left the paths differing from one
another even when handed identical innovations. All flavors now agree
across identical paths, which is checked directly.

# tsgarch 1.0.4

* Now returning the series name of the data in the spec object for use in
multivariate models.
* If the hessian cannot be calculated in the first or second step for the scaling, then
a numerical approximation is used instead. There seems to be some instability
in the egarch-nig model for some datasets.
* Added R (>= 4.1.0) requirement for use of |> pipe.
* Added additional documentation on expected index for y which must be Date or POSIXct
(not yearmon or yearqtr).
* Replaced Rf_error with (Rf_error) per Rcpp team instructions.
* Fixed the macOS C++20/TMB build issue by redefining INFINITY as a double 
infinity before including TMB.hpp in src/TMB/tsgarch_TMBExports.cpp.
* Made test fixtures lazy in tests/testthat/helper-global.R with delayedAssign(), 
so expensive model estimates are only computed when needed.
* Reduced the heavy vignette simulation workload in vignettes/demonstration.Rmd 
from 500 x 10000 to 100 x 1000.
* Cleaned 1:nrow(...) usage in edited files by switching to seq_len(nrow(...)).






# tsgarch 1.0.3

* Added the log-likelihood vector to the returned fitted and filtered object
as this will be needed in the calculation of the standard errors in the upcoming
multivariate GARCH package.
* Added an extra option to the estimation method which adds the TMB object
to the returned estimation object. This can then be used to directly vary
the parameters and quickly extract information from the filtration. This
could be also achieved by the tsfilter method with a specification object
but is much slower. Use case is for the multivariate GARCH partitioned 
hessian calculation.
* Added plus overload method to combine together GARCH specifications to generate a
multi-specification object which can then be estimated in parallel. This
is required for 2-stage multivariate GARCH models. A separate to_multi_estimate
function is also added to instead convert a list of estimated objects to
a validated multi_estimate class. Extractors include fitted, residuals and sigma.
* Removed RcppArmadillo dependency and converted code to RcppEigen since it 
is already in use by TMB.
* Switched to using simulate for the parametric simulation for the predict method.
The bootstrap still remains the most valid approach as the out of sample distribution
is best approximated by the bootstrapped residuals rather than the imposed
parametric distribution with estimated parameters.
* For the bootstrap simulated prediction, the re-sampled standardized innovations
are now scaled to avoid bias.
* Fix to h=1 and nsim. Previously when h=1, nsim was set to zero.

# tsgarch 1.0.2

* Moved a unit tests back to original folder and added a tolerance
to the expectation per CRAN maintainers directions.

# tsgarch 1.0.1

* Moved a couple of unit tests to other folder to avoid checking on CRAN
since the M1 mac had different rounding errors than other architectures
for simulation tests.


# tsgarch 1.0.0

* Initial CRAN submission.
* Changes to initialization of recursion in the ARCH equation to be more consistent
with the literature. This leads to a more than doubling in the accuracy against the
FCP benchmark.
* Correction to EGARCH forecast to account for log bias.
* Added a couple more data series for benchmarking.
* Added a pdf vignette.
* Added demo html vignettes.
* Added extensive unit tests.

