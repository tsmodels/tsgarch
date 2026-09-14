# Durbin-Levinson (Jones, 1980) reparameterization helpers for the ARMA mean
# equation. A sequence of raw coefficients r_1,...,r_n each confined to the
# open interval (-1,1) maps bijectively, via the recursion below, onto the
# stationary region of AR(n) coefficient space (Barndorff-Nielsen and Schou,
# 1973; Jones, 1980). The same mapping applied to (a) the AR coefficients and
# (b) the negated MA coefficients is used by \code{stats::arima} (via
# \sQuote{transform.pars = TRUE}) to guarantee stationarity of the AR
# polynomial and invertibility of the MA polynomial without any additional
# nonlinear inequality constraints. The forward map (pacf -> coefficients) is
# implemented identically (as a simple recursive loop) in the TMB C++
# template so that it differentiates exactly under automatic differentiation
# (see src/TMB/durbinlevinson.h).

#' Durbin-Levinson forward transform: partial autocorrelations to coefficients
#'
#' @description Maps a vector of raw parameters (each expected to lie in
#' \sQuote{(-1,1)}) to AR (or MA) coefficients that are guaranteed to lie in
#' the stationarity (invertibility) region.
#' @param r a numeric vector with values in \sQuote{(-1,1)}.
#' @return a numeric vector of the same length with the implied AR/MA
#' coefficients.
#' @keywords internal
#' @noRd
pacf_to_coef <- function(r)
{
    n <- length(r)
    if (n == 0) return(numeric(0))
    phi <- r[1]
    if (n == 1) return(phi)
    for (k in 2:n) {
        phi_new <- numeric(k)
        phi_new[k] <- r[k]
        if (k > 1) {
            phi_new[1:(k - 1)] <- phi[1:(k - 1)] - r[k] * rev(phi[1:(k - 1)])
        }
        phi <- phi_new
    }
    return(phi)
}

#' Durbin-Levinson inverse transform: coefficients to partial autocorrelations
#'
#' @description Inverse of \code{pacf_to_coef}. Used only to generate warm-start
#' values for the optimizer (e.g. from an initial \code{stats::arima} fit); the
#' forward map is the one that matters for guaranteeing stationarity/invertibility
#' during estimation.
#' @param phi a numeric vector of AR (or MA) coefficients assumed to already lie
#' in the stationarity (invertibility) region.
#' @return a numeric vector of the same length with values in \sQuote{(-1,1)}.
#' @keywords internal
#' @noRd
coef_to_pacf <- function(phi)
{
    n <- length(phi)
    if (n == 0) return(numeric(0))
    r <- numeric(n)
    r[n] <- phi[n]
    if (n > 1) {
        for (k in n:2) {
            denom <- 1 - phi[k]^2
            if (!is.finite(denom) || abs(denom) < 1e-8) {
                # fall back to a small value inside the valid region rather than
                # failing on a (numerically) non-stationary/non-invertible input
                phi <- rep(0, k - 1)
                r[k - 1] <- 0
                next
            }
            phi_new <- (phi[1:(k - 1)] + phi[k] * rev(phi[1:(k - 1)])) / denom
            r[k - 1] <- phi_new[k - 1]
            phi <- phi_new
        }
    }
    # clip for numerical safety
    pmin(pmax(r, -0.995), 0.995)
}

#' AR coefficients from raw (pacf-space) parameters
#' @keywords internal
#' @noRd
pacf_to_ar <- function(r) pacf_to_coef(r)

#' MA coefficients from raw (pacf-space) parameters
#'
#' @details Uses the AR/MA duality of the Durbin-Levinson recursion: applying
#' the transform to the negated raw parameters and negating the result yields
#' MA coefficients (in the \sQuote{y_t = e_t + theta_1 e_{t-1} + ...}
#' convention) whose characteristic polynomial has all roots outside the unit
#' circle, i.e. an invertible MA polynomial.
#' @keywords internal
#' @noRd
pacf_to_ma <- function(r) -pacf_to_coef(-r)

#' Raw (pacf-space) parameters from AR coefficients
#' @keywords internal
#' @noRd
ar_to_pacf <- function(phi) coef_to_pacf(phi)

#' Raw (pacf-space) parameters from MA coefficients
#' @keywords internal
#' @noRd
ma_to_pacf <- function(theta) -coef_to_pacf(-theta)

#' Whether an arpacf/mapacf parmatrix group is fully fixed, fully free, or absent
#'
#' @description To let users fix an entire AR or MA polynomial at target
#' coefficients (see \code{\link{garch_modelspec}}) simply by setting
#' \sQuote{value} (to the desired ar/ma coefficients, not the internal pacf
#' parameterization) and \sQuote{estimate = 0} on every row of the relevant
#' group, the two states \dQuote{fully fixed} (every row \sQuote{estimate ==
#' 0}) and \dQuote{fully free} (every row \sQuote{estimate == 1}) are given
#' different semantics for \sQuote{value} elsewhere in the package (see
#' \code{\link{arma_coefficients}} and \code{.tmb_initialize_model()}).
#' Fixing only \emph{some} of the lags in a group is not supported - the
#' Durbin-Levinson reparameterization couples all lags of a given polynomial
#' together, so there is no way to hold one lag's coefficient at a fixed
#' target while leaving others free without either giving up the
#' stationarity/invertibility guarantee or reintroducing a nonlinear
#' constraint - and is rejected with an informative error.
#' @param parmatrix a model's parmatrix \code{data.table}.
#' @param group_name either \sQuote{"arpacf"} or \sQuote{"mapacf"}.
#' @return one of \sQuote{"none"} (the group has no rows, i.e. that
#' polynomial's order is zero), \sQuote{"fixed"} (every row has
#' \sQuote{estimate == 0}) or \sQuote{"free"} (every row has
#' \sQuote{estimate == 1}).
#' @keywords internal
#' @noRd
arma_block_status <- function(parmatrix, group_name)
{
    group <- estimate <- NULL
    rows <- parmatrix[group == group_name]
    if (NROW(rows) == 0) return("none")
    est <- rows$estimate
    if (length(unique(est)) > 1) {
        lag_label <- if (group_name == "arpacf") "ar" else "ma"
        stop("\narma: partial fixing of individual ", lag_label,
             " lags is not supported (the Durbin-Levinson reparameterization\n",
             "couples all lags of a given polynomial together). Either fix all ",
             group_name, " parameters (estimate = 0, value = the desired ",
             lag_label, " coefficients) or leave them all free (estimate = 1).")
    }
    if (all(est == 0)) "fixed" else "free"
}

#' Validate a fixed target AR/MA polynomial is stationary/invertible
#'
#' @description Called whenever an entire AR or MA polynomial is fixed (see
#' \code{\link{arma_block_status}}) to check, before transforming to the
#' internal pacf parameterization via \code{\link{ar_to_pacf}}/
#' \code{\link{ma_to_pacf}}, that the user-supplied target coefficients
#' actually correspond to a stationary (AR) or invertible (MA) polynomial;
#' \code{coef_to_pacf}'s internal numerical fallback for near-degenerate
#' inputs would otherwise silently substitute a different value rather than
#' erroring on a genuinely invalid target.
#' @param target numeric vector of target AR (or MA) coefficients.
#' @param type either \sQuote{"ar"} or \sQuote{"ma"}.
#' @keywords internal
#' @noRd
validate_arma_fixed_target <- function(target, type = c("ar","ma"))
{
    type <- match.arg(type)
    if (length(target) == 0) return(invisible(TRUE))
    if (type == "ar") {
        roots <- polyroot(c(1, -target))
        label <- "AR"; requirement <- "stationary"
    } else {
        roots <- polyroot(c(1, target))
        label <- "MA"; requirement <- "invertible"
    }
    if (!all(is.finite(roots)) || min(Mod(roots)) <= 1) {
        article <- if (requirement == "invertible") "an" else "a"
        stop("\narma: the fixed ", label, " coefficients (", paste(signif(target, 4), collapse = ", "),
             ") do not correspond to ", article, " ", requirement, " polynomial (a root lies within or on the unit circle).")
    }
    invisible(TRUE)
}

#' Construct the full [M] parmatrix block for a model's mean equation
#'
#' @description Builds the \sQuote{arpacf}/\sQuote{mapacf}/\sQuote{tau}
#' parameter rows shared by every \sQuote{.parameters_<model>} initializer
#' (see R/initialization.R). The arpacf/mapacf rows hold the raw
#' Durbin-Levinson (partial autocorrelation) parameters rather than the AR/MA
#' coefficients themselves; any value in \sQuote{(-1,1)} maps (inside the TMB
#' template) to AR coefficients in the stationarity region and MA coefficients
#' in the invertibility region, so no additional nonlinear constraint is
#' required. The actual ar/ma coefficients implied by these raw values are
#' available separately via \code{\link{arma_coefficients}}. The
#' \sQuote{tau} rows hold the mean-equation regressor coefficients
#' (\sQuote{tau1..taum} when \code{xreg} is present, one fixed dummy row
#' otherwise, matching the arpacf/mapacf dummy convention and keeping the
#' parmatrix row count aligned with the always >= 1 column count of the
#' regressor matrix passed to TMB). All starting values (ar/ma pacf and tau)
#' come from a single joint \code{stats::arima(..., xreg = )} fit so they are
#' mutually consistent. The +/-100 bounds on tau match the existing xi
#' handling; users whose regressors are on a very different scale from y may
#' need to widen these bounds in the spec's parmatrix.
#' @param y numeric vector, the (not necessarily demeaned) target series;
#' only used, together with \code{mu}, to generate \code{stats::arima}-based
#' starting values.
#' @param mu the (scalar) constant/unconditional mean value.
#' @param arma length 2 integer vector \sQuote{c(ar, ma)}.
#' @param xreg optional matrix of mean-equation regressors (ncol = number of
#' \sQuote{tau} parameters); \code{NULL} for none.
#' @return a \code{data.table} with the same columns as the rest of a
#' \sQuote{parmatrix} (parameter, value, lower, upper, estimate, scale,
#' group, equation, symbol), with \sQuote{ar_order} (or 1, as a fixed dummy
#' row when \sQuote{ar_order = 0}) rows of group \sQuote{arpacf}, followed by
#' \sQuote{ma_order} (or 1 dummy) rows of group \sQuote{mapacf}, followed by
#' \sQuote{ncol(xreg)} (or 1 dummy) rows of group \sQuote{tau}.
#' @keywords internal
#' @noRd
arma_parmatrix_rows <- function(y, mu, arma, xreg = NULL)
{
    ar_order <- arma[1]
    ma_order <- arma[2]
    if (sum(arma) > 0 || !is.null(xreg)) {
        arma_pacf <- initialize_arma_pacf(y - mu, arma, xreg)
    } else {
        arma_pacf <- list(ar = numeric(0), ma = numeric(0), tau = NULL)
    }
    if (ar_order == 0) {
        ar_rows <- data.table("parameter" = "arpacf1", value = 0,
                              lower = -0.995, upper = 0.995, estimate = 0,
                              scale = 1, group = "arpacf", equation = "[M]",
                              symbol = "r^{ar}_1")
    } else {
        ar_rows <- data.table("parameter" = paste0("arpacf",1:ar_order),
                              value = arma_pacf$ar, lower = -0.995, upper = 0.995,
                              estimate = 1, scale = 1, group = "arpacf",
                              equation = "[M]",
                              symbol = paste0("r^{ar}_",1:ar_order))
    }
    if (ma_order == 0) {
        ma_rows <- data.table("parameter" = "mapacf1", value = 0,
                              lower = -0.995, upper = 0.995, estimate = 0,
                              scale = 1, group = "mapacf", equation = "[M]",
                              symbol = "r^{ma}_1")
    } else {
        ma_rows <- data.table("parameter" = paste0("mapacf",1:ma_order),
                              value = arma_pacf$ma, lower = -0.995, upper = 0.995,
                              estimate = 1, scale = 1, group = "mapacf",
                              equation = "[M]",
                              symbol = paste0("r^{ma}_",1:ma_order))
    }
    if (is.null(xreg)) {
        tau_rows <- data.table("parameter" = "tau1", value = 0,
                               lower = -100, upper = 100, estimate = 0,
                               scale = 1, group = "tau", equation = "[M]",
                               symbol = "\\tau_1")
    } else {
        m <- NCOL(xreg)
        tau_rows <- data.table("parameter" = paste0("tau",1:m),
                               value = arma_pacf$tau, lower = -100, upper = 100,
                               estimate = 1, scale = 1, group = "tau",
                               equation = "[M]",
                               symbol = paste0("\\tau_",1:m))
    }
    rbind(ar_rows, ma_rows, tau_rows)
}

#' Initial values for the ARMA mean equation parameters
#'
#' @description Fits an unconstrained \code{stats::arima} model to (demeaned)
#' \code{y} to obtain reasonable starting values for the raw (pacf-space) ar/ma
#' parameters used internally by the TMB Durbin-Levinson mean recursion.
#' @param y numeric vector (already demeaned if a constant is separately
#' estimated).
#' @param arma length 2 integer vector \sQuote{c(ar, ma)}.
#' @param xreg optional matrix of mean-equation regressors; when non-NULL a
#' single joint \code{stats::arima(y, xreg = xreg)} fit supplies mutually
#' consistent warm starts for both the ar/ma pacf block and the tau block.
#' @return a list with elements \sQuote{ar}, \sQuote{ma} and \sQuote{tau}:
#' the first two numeric vectors of raw (pacf-space) starting values bounded
#' away from +/-1, \sQuote{tau} a numeric vector of regressor-coefficient
#' starting values (\code{NULL} when \code{xreg} is \code{NULL}).
#' @keywords internal
#' @noRd
initialize_arma_pacf <- function(y, arma, xreg = NULL)
{
    ar <- arma[1]
    ma <- arma[2]
    ar_pacf <- if (ar > 0) rep(0.01, ar) else numeric(0)
    ma_pacf <- if (ma > 0) rep(0.01, ma) else numeric(0)
    tau <- NULL
    has_xreg <- !is.null(xreg)
    if (has_xreg) {
        xreg <- coredata(xreg)
        # OLS fallback for the tau starts (also the source of the tau starts
        # when sum(arma) == 0 and no arima fit is run at all)
        ols <- try(qr.solve(crossprod(xreg), crossprod(xreg, as.numeric(y))), silent = TRUE)
        tau <- if (!inherits(ols, "try-error") && all(is.finite(ols))) as.numeric(ols) else rep(0, NCOL(xreg))
    }
    if (sum(arma) == 0) {
        return(list(ar = ar_pacf, ma = ma_pacf, tau = tau))
    }
    fit <- try(stats::arima(as.numeric(y), order = c(ar, 0, ma), include.mean = FALSE,
                            xreg = if (has_xreg) xreg else NULL,
                            method = "ML", transform.pars = TRUE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
        cf <- stats::coef(fit)
        if (ar > 0) {
            ar_coef <- as.numeric(cf[paste0("ar", 1:ar)])
            if (all(is.finite(ar_coef))) {
                tmp <- try(ar_to_pacf(ar_coef), silent = TRUE)
                if (!inherits(tmp, "try-error") && all(is.finite(tmp))) ar_pacf <- tmp
            }
        }
        if (ma > 0) {
            ma_coef <- as.numeric(cf[paste0("ma", 1:ma)])
            if (all(is.finite(ma_coef))) {
                tmp <- try(ma_to_pacf(ma_coef), silent = TRUE)
                if (!inherits(tmp, "try-error") && all(is.finite(tmp))) ma_pacf <- tmp
            }
        }
        if (has_xreg) {
            # arima names the regressor coefficients after the xreg colnames,
            # which check_xreg() guarantees are set (x1..xm when absent);
            # index them positionally (the trailing ar+ma coefficients)
            # rather than by name for robustness. For xreg_type = "armax"
            # this is the arma_errors-form estimate and serves only as a
            # warm start - no conversion is attempted.
            tau_cf <- as.numeric(tail(cf, NCOL(xreg)))
            if (all(is.finite(tau_cf))) tau <- tau_cf
        }
    }
    list(ar = ar_pacf, ma = ma_pacf, tau = tau)
}

#' Extend the ARMA conditional mean with newly observed data
#'
#' @description Continues the ARMA mean recursion
#' \sQuote{conditional_mean_t = mu + sum_i phi_i*(y_{t-i}-mu) + sum_j theta_j*eps_{t-j}}
#' (the same recursion used inside the TMB template, see
#' \code{src/TMB/garchfun.hpp}, where \sQuote{eps_t = y_t - conditional_mean_t})
#' for newly appended observations, using the already-computed historical
#' conditional mean as the continuation state (residuals are always
#' recoverable as \sQuote{y - conditional_mean}, see
#' \code{fitted}/\code{residuals} methods). This lets \code{\link{tsfilter}}
#' incrementally update the fitted mean for new data without re-running the
#' full TMB likelihood.
#' @param y_full numeric vector of the full (old + new), merged series.
#' @param n_old integer, the number of observations already covered by
#' \code{old_conditional_mu} (i.e. \sQuote{length(y_full) - n_old} new
#' observations will be computed).
#' @param mu the (scalar) unconditional mean.
#' @param ar numeric vector of AR coefficients (already transformed via
#' \code{\link{pacf_to_ar}}; may be length zero).
#' @param ma numeric vector of MA coefficients (already transformed via
#' \code{\link{pacf_to_ma}}; may be length zero).
#' @param old_conditional_mu numeric vector of length \sQuote{n_old} with the
#' already-computed historical conditional mean.
#' @return a numeric vector of length \sQuote{length(y_full) - n_old} with the
#' conditional mean for the newly appended observations only.
#' @keywords internal
#' @noRd
arma_filter_extend <- function(y_full, n_old, mu, ar, ma, old_conditional_mu)
{
    ar_order <- length(ar)
    ma_order <- length(ma)
    n_new <- length(y_full) - n_old
    if (n_new <= 0) return(numeric(0))
    y_full <- as.numeric(y_full)
    z <- y_full - mu
    old_eps <- head(y_full, n_old) - old_conditional_mu
    eps <- c(old_eps, rep(0, n_new))
    cond_mu <- c(old_conditional_mu, rep(mu, n_new))
    for (i in seq_len(n_new)) {
        t <- n_old + i
        mean_t <- 0
        if (ar_order > 0) {
            for (j in seq_len(ar_order)) {
                idx <- t - j
                mean_t <- mean_t + ar[j] * (if (idx <= 0) 0 else z[idx])
            }
        }
        if (ma_order > 0) {
            for (j in seq_len(ma_order)) {
                idx <- t - j
                mean_t <- mean_t + ma[j] * (if (idx <= 0) 0 else eps[idx])
            }
        }
        eps[t] <- z[t] - mean_t
        cond_mu[t] <- mu + mean_t
    }
    cond_mu[(n_old + 1):(n_old + n_new)]
}

#' Inverse AR/MA Roots
#'
#' @description Computes the inverse roots of the AR and MA polynomials of the
#' mean equation in a fitted or specified ARMA-GARCH model.
#' @param object an object of class \dQuote{tsgarch.estimate} or
#' \dQuote{tsgarch.spec}.
#' @return A list with elements \sQuote{ar} and \sQuote{ma}, each a complex
#' vector of inverse roots. The AR polynomial is
#' \eqn{1 - \phi_1 z - \dots - \phi_p z^p}{1 - phi_1 z - ... - phi_p z^p} and
#' the MA polynomial is
#' \eqn{1 + \theta_1 z + \dots + \theta_q z^q}{1 + theta_1 z + ... + theta_q z^q}.
#' If either component is absent, the corresponding vector is \code{complex(0)}.
#' @export
arma_inverse_roots <- function(object)
{
    ac <- arma_coefficients(object)
    out <- list(ar = complex(0), ma = complex(0))
    if (length(ac$ar) > 0) {
        arpoly <- c(1, -as.numeric(ac$ar))
        r <- polyroot(arpoly)
        out$ar <- 1/r
        out$ar[Mod(r) == 0] <- NaN + 1i * NaN
    }
    if (length(ac$ma) > 0) {
        mapoly <- c(1, as.numeric(ac$ma))
        r <- polyroot(mapoly)
        out$ma <- 1/r
        out$ma[Mod(r) == 0] <- NaN + 1i * NaN
    }
    return(out)
}

#' ARMA Impulse Response Function
#'
#' @description Computes the impulse response function (IRF) for the mean
#' equation of a fitted or specified ARMA-GARCH model.
#' @param object an object of class \dQuote{tsgarch.estimate} or
#' \dQuote{tsgarch.spec}.
#' @param lag.max integer, the maximum lag. If \code{NULL}, an adaptive
#' horizon is chosen by truncating when the response falls below
#' \code{1e-3} of the peak absolute value, with a floor of 10 and a cap of 100.
#' @return A list with elements \sQuote{psi} (the IRF, indexed from lag 0) and
#' \sQuote{cumulative} (the cumulative sum).
#' @export
arma_irf <- function(object, lag.max = NULL)
{
    ac <- arma_coefficients(object)
    ar <- as.numeric(ac$ar)
    ma <- as.numeric(ac$ma)
    L <- lag.max
    if (is.null(L)) {
        L <- 200
    }
    psi <- c(1, as.numeric(stats::ARMAtoMA(ar = ar, ma = ma, lag.max = L)))
    if (is.null(lag.max)) {
        threshold <- 1e-3 * max(abs(psi))
        last <- max(which(abs(psi) >= threshold))
        last <- min(max(last, 10), 100)
        psi <- psi[1:(last + 1)]
    }
    cum <- cumsum(psi)
    names(psi) <- paste0("lag", seq_along(psi) - 1)
    return(list(psi = psi, cumulative = cum))
}

#' Near AR/MA Root Cancellation
#'
#' @description Identifies AR and MA inverse roots that are close in the
#' complex plane, a visual cue for near-common-factors in the mean equation.
#' @param roots a list, normally the output of \code{\link{arma_inverse_roots}}.
#' @param tol numeric distance threshold in the complex plane.
#' @return A \code{data.table} with one row per near-cancelling pair and columns
#' \sQuote{ar_index}, \sQuote{ma_index}, \sQuote{ar_root}, \sQuote{ma_root}
#' and \sQuote{distance}.
#' @export
arma_near_cancellation <- function(roots, tol = 0.1)
{
    ar <- as.complex(roots$ar)
    ma <- as.complex(roots$ma)
    out <- data.table::data.table(ar_index = integer(0), ma_index = integer(0),
                                  ar_root = complex(0), ma_root = complex(0),
                                  distance = numeric(0))
    if (length(ar) == 0 || length(ma) == 0) return(out)
    dm <- sapply(ar, function(a) Mod(a - ma))
    if (is.matrix(dm)) {
        idx <- which(dm < tol, arr.ind = TRUE)
        if (nrow(idx) > 0) {
            for (k in seq_len(nrow(idx))) {
                i <- idx[k, 1]; j <- idx[k, 2]
                out <- rbind(out, data.table::data.table(
                    ar_index = i, ma_index = j,
                    ar_root = ar[i], ma_root = ma[j],
                    distance = dm[i, j]), use.names = TRUE)
            }
        }
    } else if (length(dm) > 0) {
        # length(ar)==1 or length(ma)==1 path, dm is a vector
        idx <- which(dm < tol)
        for (k in idx) {
            i <- ifelse(length(ar) > 1, k, 1)
            j <- ifelse(length(ma) > 1, k, 1)
            out <- rbind(out, data.table::data.table(
                ar_index = i, ma_index = j,
                ar_root = ar[i], ma_root = ma[j],
                distance = dm[k]), use.names = TRUE)
        }
    }
    return(out)
}

#' Parametric-uncertainty draws of the standardized residuals
#'
#' @description Draws \code{B} perturbed parameter vectors from the
#' asymptotic normal approximation \eqn{\hat\theta + N(0, V)}, where \eqn{V}
#' is \code{\link{vcov}}, and for each draw re-filters the \emph{actual
#' observed data} (no re-simulation, no re-optimization) via a cheap forward
#' TMB \code{report()} call to obtain the implied conditional mean and
#' conditional standard deviation, and hence a perturbed standardized
#' residual series \eqn{z_t = (y_t - \mu_{t,b})/\sigma_{t,b}}. This
#' characterizes how much the standardized-residual diagnostics themselves
#' wobble due to parameter estimation uncertainty. Perturbed draws are
#' clipped to each parameter's \code{lower}/\code{upper} bounds (as stored in
#' \code{object$parmatrix}) to keep the forward pass numerically valid; this
#' is a deliberate simplification appropriate for a diagnostic plot, not a
#' formal inferential procedure.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param B integer, the number of draws.
#' @param vcov_type the type of covariance matrix to use (see \code{\link{vcov}}).
#' @return A numeric matrix with \code{B} columns, one perturbed standardized
#' residual series per column, and \code{length(residuals(object))} rows.
#' @keywords internal
#' @noRd
.parametric_standardized_residual_draws <- function(object, B = 500, vcov_type = "H")
{
    estimate <- NULL
    theta_hat <- coef(object)
    k <- length(theta_hat)
    V <- vcov(object, type = vcov_type)
    Lchol <- tryCatch(chol(V), error = function(e) {
        tryCatch(chol(V + diag(1e-8, k)),
                 error = function(e2) stop("\nvcov(object, type = '", vcov_type, "') is not positive-definite; try a different vcov_type or envelope = 'simulate'."))
    })
    lower <- object$parmatrix[estimate == 1]$lower
    upper <- object$parmatrix[estimate == 1]$upper

    full_spec <- object$spec
    full_spec$parmatrix <- data.table::copy(object$parmatrix)
    model_init <- .tmb_initialize_model(full_spec)
    tmb <- TMB::MakeADFun(data = model_init$data, parameters = model_init$parameters,
                          map = model_init$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    m <- full_spec$model_options[1]
    y <- as.numeric(object$spec$target$y_orig)
    n <- length(y)

    z <- matrix(NA_real_, nrow = n, ncol = B)
    for (b in seq_len(B)) {
        theta_b <- as.numeric(theta_hat) + as.numeric(t(Lchol) %*% stats::rnorm(k))
        theta_b <- pmin(pmax(theta_b, lower), upper)
        rep_b <- tmb$report(theta_b)
        sigma_b <- rep_b$sigma
        mu_b <- rep_b$conditional_mean
        if (m > 0) {
            sigma_b <- sigma_b[-seq_len(m)]
            mu_b <- mu_b[-seq_len(m)]
        }
        z[, b] <- (y - mu_b) / sigma_b
    }
    return(z)
}
