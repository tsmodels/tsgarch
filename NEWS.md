# tsgarch 1.0.5

* Added a jointly estimated ARMA(p,q) mean equation, via a new `arma` argument
to `garch_modelspec` (defaults to `c(0,0)`, fully backward compatible).
Available for all 8 GARCH flavors, including `igarch` and `ewma`. The
`constant` argument remains independent of `arma` and continues to control
only whether the unconditional mean `mu` is estimated or fixed at zero.
* The AR/MA coefficients are guaranteed stationary/invertible by construction,
via a Durbin-Levinson/Jones partial-autocorrelation (PACF) reparameterization
of the raw optimization parameters, so no additional nonlinear constraints
are required; exact (autodiff, not finite-difference) Jacobians of the ARMA
mean equation are available throughout via TMB. The transformed AR/MA
coefficients can be extracted with the new `arma_coefficients()` function.
* `fitted()` now returns the genuinely time-varying conditional mean when an
ARMA mean equation is specified (previously always a constant `mu`, now
completing the documented "vector the size of y" contract), and `residuals()`
is now always exactly `y - fitted(object)` for every model, including a fix
for models with no ARMA (a small, purely internal simplification with no
behavior change there).
* `predict()`, `simulate()` and `tsfilter()` are all ARMA-aware: point
forecasts, simulated paths and incremental filtering all correctly account
for the AR/MA dynamics in the mean equation. The combined pre-sample/burn-in
length used internally is `max(garch order, arma order)`.
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

