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

#' Initial values for the ARMA mean equation parameters
#'
#' @description Fits an unconstrained \code{stats::arima} model to (demeaned)
#' \code{y} to obtain reasonable starting values for the raw (pacf-space) ar/ma
#' parameters used internally by the TMB Durbin-Levinson mean recursion.
#' @param y numeric vector (already demeaned if a constant is separately
#' estimated).
#' @param arma length 2 integer vector \sQuote{c(ar, ma)}.
#' @return a list with elements \sQuote{ar} and \sQuote{ma}, each a numeric
#' vector of raw (pacf-space) starting values bounded away from +/-1.
#' @keywords internal
#' @noRd
initialize_arma_pacf <- function(y, arma)
{
    ar <- arma[1]
    ma <- arma[2]
    ar_pacf <- if (ar > 0) rep(0.01, ar) else numeric(0)
    ma_pacf <- if (ma > 0) rep(0.01, ma) else numeric(0)
    if (sum(arma) == 0) {
        return(list(ar = ar_pacf, ma = ma_pacf))
    }
    fit <- try(stats::arima(as.numeric(y), order = c(ar, 0, ma), include.mean = FALSE,
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
    }
    list(ar = ar_pacf, ma = ma_pacf)
}

#' Extend the ARMA mean equation residuals with newly observed data
#'
#' @description Continues the ARMA mean recursion
#' \sQuote{(y_t - mu) = sum_i phi_i*(y_{t-i}-mu) + eps_t + sum_j theta_j*eps_{t-j}}
#' (the same recursion used inside the TMB template, see
#' \code{src/TMB/garchfun.hpp}) for newly appended observations, using the
#' already-computed historical residuals as the continuation state. This
#' lets \code{\link{tsfilter}} incrementally update residuals for new data
#' without re-running the full TMB likelihood.
#' @param y_full numeric vector of the full (old + new), merged series.
#' @param n_old integer, the number of observations already covered by
#' \code{old_residuals} (i.e. \sQuote{length(y_full) - n_old} new observations
#' will be computed).
#' @param mu the (scalar) unconditional mean.
#' @param ar numeric vector of AR coefficients (already transformed via
#' \code{\link{pacf_to_ar}}; may be length zero).
#' @param ma numeric vector of MA coefficients (already transformed via
#' \code{\link{pacf_to_ma}}; may be length zero).
#' @param old_residuals numeric vector of length \sQuote{n_old} with the
#' already-computed historical residuals.
#' @return a numeric vector of length \sQuote{length(y_full) - n_old} with the
#' residuals for the newly appended observations only.
#' @keywords internal
#' @noRd
arma_filter_extend <- function(y_full, n_old, mu, ar, ma, old_residuals)
{
    ar_order <- length(ar)
    ma_order <- length(ma)
    n_new <- length(y_full) - n_old
    if (n_new <= 0) return(numeric(0))
    z <- as.numeric(y_full) - mu
    eps <- c(old_residuals, rep(0, n_new))
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
    }
    eps[(n_old + 1):(n_old + n_new)]
}
