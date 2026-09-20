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
* `predict()`, `simulate()`, `tsfilter()`, `garch_backtest()` and
`tsprofile()` all account for the full ARMA-X mean equation: point
forecasts, both the parametric and bootstrap predictive distributions,
simulated paths, incremental filtering, each rolling window's refit and the
simulate/re-estimate profile all reflect the AR/MA dynamics and the
regressor contribution. Variance targeting uses the unconditional variance
of the ARMA innovations rather than deviations from the constant `mu`, so
the R side `unconditional()`/`target_omega` agree with the TMB side. Mean
regressors are refused by `tsprofile()`, as `vreg` already was, since its
simulation step does not carry them into the simulated sample. The combined
pre-sample/burn-in length used internally is `max(garch order, arma order)`.
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
* `tsequation()`, and the `as_flextable` summary footer built from it, now
render the conditional mean equation `eq_mean` first in the equation block,
covering the constant, the ARMA terms and the mean regressors under both
`xreg_type` conventions.
* Fixed `tsbacktest()` silently refitting the wrong model inside its rolling
windows: the inner `garch_modelspec()` call passed no `model`, so e.g. an
`egarch` spec was backtested as a vanilla `garch`. The spec now records the
user-facing model name in `model$model_name`, necessary because `ewma` is
coerced to `igarch` internally, and each window is refit with it.
* Fixed spec re-specification dropping the `ewma` restriction. `ewma` is
coerced to `igarch` internally, the difference being only the fixed `omega`
row of the `parmatrix`, and both `.spec2newspec()` and `garch_profile()`
rebuilt specs from the coerced name. A filtered `ewma` model therefore
carried a spurious free `omega`, leaving `npars` and AIC/BIC out by one
parameter's worth, and `garch_profile()` re-fit an `igarch` model. Filtered
`sigma` is unchanged; filtered `ewma` AIC/BIC shift accordingly.
* `simulate()` now raises an error instead of silently returning zeros when
the implied initial variance is not positive and finite. For an `ewma`
specification `omega` is fixed at zero, so the seed was zero and the
variance recursion stayed there at every step, returning `sigma` identically
zero and a series equal to `mu`. Supply `var_init` for such models;
`predict()` always did, and every other flavor has a positive implied seed
and is unaffected.
* `simulate()` on a spec taken from an estimated object (`mod$spec`, which
carries `parmatrix = NULL`) no longer aborts the R session, raising a clean
error asking for the estimated `parmatrix` to be assigned onto the spec
first. A specification serialized by a version up to 1.0.4 predates the
`arpacf`, `mapacf` and `tau` rows of the `parmatrix` and is likewise refused
with a message naming `garch_modelspec()` as the remedy, rather than failing
with `Error when reading the variable: 'arpacf'`.
* Fixed `tsfilter()`'s zero fill for a missing `newvreg` being built with
`as.matrix(0, ...)`, a 1x1 matrix, instead of `matrix(0, nrow, ncol)`.
* Corrected the pre-sample initialization of the ARCH recursion for the
asymmetric flavors (`egarch`, `gjrgarch`, `aparch`, `fgarch`). Whether a lag
lookback falls inside the pre-sample was tested against the ARCH order
rather than against the pre-sample length, and the two coincide only when
the ARCH order is the larger. Observations just past the pre-sample
therefore took a zeroed residual in place of the initial ARCH value, in the
likelihood and in `simulate()` alike; both now test the pre-sample length.
This is reachable through a GARCH order with `q > p`, so `logLik` and the
coefficients shift for these four flavors in that case. `order = c(1,1)` is
arithmetically unchanged, and `garch`, `igarch` and `ewma` are unaffected at
every order, their pre-sample ARCH input being constant.
* Corrected the pre-sample initialization of `simulate()`, which now
reproduces a fitted model's `sigma` to machine precision for every flavor,
where the discrepancy previously reached 5e-2. `var_init` and `innov_init`
are documented as seeding every sample path identically, but were filled
column-major for several flavors, rotating them across the paths instead of
broadcasting them; the lag-indexed ARCH initialization was paired with the
oldest pre-sample period rather than the most recent, mismatching each lag
with its coefficient; `aparch` and `fgarch` recycled their per-lag `gamma`
(and `fgarch` its `eta`) across paths rather than across lags; and an
`arch_initial` shorter than the pre-sample was collapsed onto its first
element rather than padded, putting the lag 1 value into every lag. These
are observable with more than one sample path, or a pre-sample longer than
one period.
* Corrected the default ARCH initialization of a simulation, which is the
expectation of each flavor's ARCH equation. `egarch` instead evaluated its
equation at the zeroed pre-sample innovations literally, contributing
`-gamma_j kappa` where both the expectation and the likelihood give zero,
which biased the opening steps downwards. `aparch` raised its already
power-transformed initial variance to `delta/2` a second time, the identity
only at `delta = 2`, so a default simulation opened away from its own
unconditional level by 3 to 10 percent for any other `delta`. A default
simulation now opens exactly at the unconditional level for every flavor
that has one.
* A pre-sample innovation supplied to an `egarch` simulation through
`innov_init` now carries its leverage effect; only its magnitude was used,
the `alpha_j z_{t-j}` term of the ARCH equation being absent from the
pre-sample. The `egarch` predictive distribution and the
simulation-approximated higher order forecast both move, and at a horizon of
one the mean of the simulated distribution now agrees with the closed form
forecast to machine precision, where it was out by around 3 percent.
* `tsbacktest()` no longer warns that its iterations drew random numbers
without declaring a seed. They do draw them, through `predict()`, so they
are now given parallel-safe streams.
* A standard error that is legitimately undefined is now reported as `NaN`
without the low-level warning behind it. When a parameter is estimated at
one of its bounds, or for `ewma`, whose unit persistence constraint leaves
its two coefficients perfectly dependent, the Hessian has no positive
definite direction there and the delta method variance is negative.
`kkt1`/`kkt2` in the estimated object's `conditions` remain the signal that
the solution sits on a bound.
* Documented the matrix form of `var_init` accepted by `simulate()` for a
`cgarch` model: a `max(order, arma)` by 2 matrix whose first column
initializes the permanent (long run) component and whose second column
initializes the total conditional variance.

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

