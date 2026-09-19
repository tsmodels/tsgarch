#' ARMA mean overlay for a simulated GARCH variance/innovation process
#'
#' @description Shared by every \sQuote{.simulate_<model>} function: applies
#' the ARMA mean recursion (\code{\link{.armasimvec}} in
#' src/simulation.cpp) on top of an already-simulated GARCH variance/
#' innovation process (\sQuote{simc$epsilon}, unmodified), mirroring
#' rugarch's own two-stage design (variance via e.g. sgarchsimC, mean
#' overlay via armaxsim using the same innovations - see rugarch's
#' src/garchsim.cpp: msgarchsim() + marmaxsim()). Returns \sQuote{simc$series}
#' unchanged when \sQuote{sum(arma_order) == 0}.
#' @param object the (spec-like) model object, as passed into each
#' \sQuote{.simulate_<model>} function.
#' @param simc the list returned by the model's Rcpp \sQuote{*simvec}
#' function, with (at least) \sQuote{series} and \sQuote{epsilon} elements.
#' @param mu the (scalar) unconditional mean.
#' @param maxpq the combined pre-sample length (\sQuote{max(garch order, arma order)}).
#' @param arma_order length 2 integer vector \sQuote{c(ar, ma)}.
#' @param nsim number of simulated paths (rows).
#' @param extra_args the \sQuote{list(...)} of extra arguments passed to the
#' calling \sQuote{.simulate_<model>} function; \sQuote{series_init}/
#' \sQuote{resid_init} (each of length \sQuote{maxpq}), if present, seed the
#' ARMA pre-sample from real history (see \code{predict.R}'s
#' \sQuote{simulated_distribution()}) rather than the unconditional mean /
#' zero shocks default.
#' @param xreg optional matrix of mean equation regressors with
#' \sQuote{h + burn} rows (already validated by the calling
#' \sQuote{.simulate_<model>} function); multiplied internally by tau.
#' @return the (possibly ARMA-adjusted) simulated series matrix.
#' @keywords internal
#' @noRd
.arma_simulate_overlay <- function(object, simc, mu, maxpq, arma_order, nsim, extra_args = list(), xreg = NULL)
{
    group <- NULL
    include_xreg <- isTRUE(object$xreg$include_xreg)
    armax <- isTRUE(object$xreg$xreg_type == "armax")
    if (sum(arma_order) == 0 && !include_xreg) {
        return(simc$series)
    }
    arma_coef <- arma_coefficients(object)
    ar <- as.numeric(arma_coef$ar)
    ma <- as.numeric(arma_coef$ma)
    series_sim <- simc$series
    # xtau covers all T = maxpq + h columns: the pre-sample columns take the
    # actual last in-sample x'tau values (the same dates as a supplied
    # series_init), the remaining columns the supplied regressor values
    # (zeros when the model has no regressors)
    tau <- object$parmatrix[group == "tau"]$value
    xreg_old <- object$xreg$xreg
    if (is.null(xreg_old)) xreg_old <- matrix(0, ncol = 1, nrow = max(1, NROW(object$target$y_orig)))
    if (length(tau) != NCOL(xreg_old)) tau <- rep(0, NCOL(xreg_old))
    xtau_insample <- if (maxpq > 0) as.numeric(tail(as.matrix(xreg_old), maxpq) %*% tau) else numeric(0)
    if (include_xreg && !is.null(xreg)) {
        xtau_new <- as.numeric(as.matrix(xreg) %*% tau)
    } else {
        xtau_new <- rep(0, ncol(series_sim) - maxpq)
    }
    xtau_full <- c(xtau_insample, xtau_new)
    # seed the pre-sample series so the AR lookback term vanishes there:
    # under arma_errors that is mu + x'tau (w = y - mu - x'tau = 0), under
    # armax it is mu (y - mu = 0) - unless the caller supplies the actual
    # last observed values via extra_args$series_init, which are used as-is
    # under both conventions.
    if (maxpq > 0) {
        if (!is.null(extra_args$series_init)) {
            if (length(extra_args$series_init) != maxpq) stop(paste0("\nseries_init must be of length max(garch order, arma order) : ", maxpq))
            series_sim[,seq_len(maxpq)] <- matrix(extra_args$series_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
        } else {
            seed <- if (armax) rep(mu, maxpq) else mu + xtau_insample
            series_sim[,seq_len(maxpq)] <- matrix(seed, ncol = maxpq, nrow = nsim, byrow = TRUE)
        }
    }
    # the pre-sample columns of simc$epsilon may default to a deterministic,
    # non-zero value (e.g. z = 1) used to seed the ARCH init term as sigma^2;
    # this is fine for a squared ARCH term but would leak a spurious
    # deterministic "shock" into the (sign-sensitive) MA lookback below, so
    # use a copy that defaults to zeroed pre-sample entries for that purpose
    # instead (unconditional start, zero shocks), unless the caller supplies
    # the actual last observed residuals via extra_args$resid_init. The
    # actual, real-valued forecast/simulated epsilon columns are unaffected
    # either way.
    ma_epsilon <- simc$epsilon
    if (maxpq > 0) {
        if (!is.null(extra_args$resid_init)) {
            if (length(extra_args$resid_init) != maxpq) stop(paste0("\nresid_init must be of length max(garch order, arma order) : ", maxpq))
            ma_epsilon[,seq_len(maxpq)] <- matrix(extra_args$resid_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
        } else {
            ma_epsilon[,seq_len(maxpq)] <- 0
        }
    }
    .armasimvec(series_sim = series_sim, epsilon = ma_epsilon, ar = ar, ma = ma, mu = mu,
                xtau = xtau_full, armax = as.integer(armax), presample = maxpq)
}

# init.col(j) in the simulation recursions is indexed by ARCH lag j, not by
# pre-sample column, so only the first order[1] entries are ever read and
# the remaining columns merely fill the pre-sample block. A short per-lag
# vector is therefore padded rather than collapsed onto its first entry,
# which would put lag 1's initialization into every lag slot. Anything
# longer than maxpq (the nsim x p matrix derived from innov_init) keeps the
# previous behaviour, since it is not a per-lag vector.
.expand_arch_initial <- function(init, maxpq, nsim)
{
    init <- as.numeric(init)
    if (length(init) > maxpq) {
        init <- rep(init[1], maxpq)
    } else if (length(init) < maxpq) {
        init <- c(init, rep(init[length(init)], maxpq - length(init)))
    }
    matrix(init, ncol = maxpq, nrow = nsim, byrow = TRUE)
}

# validate the simulate() mean-regressor argument: a matrix/xts of h + burn
# rows by NCOL(xreg) columns (deliberately unlike the legacy vreg argument,
# which is a pre-multiplied vector - see AGENTS.md). NULL with regressors in
# the model warns and substitutes zeros.
.process_simulation_xreg <- function(object, xreg, h)
{
    include_xreg <- isTRUE(object$xreg$include_xreg)
    ncols <- if (!is.null(object$xreg$xreg)) NCOL(object$xreg$xreg) else 1L
    if (include_xreg) {
        if (is.null(xreg)) {
            warning("\nxreg is NULL but the model was specified with mean regressors; setting to zero.")
            xreg <- matrix(0, nrow = h, ncol = ncols)
        }
        if (!is.xts(xreg)) xreg <- as.matrix(xreg)
        if (NROW(xreg) != h) stop("\nxreg must have h + burn rows.")
        if (NCOL(xreg) != ncols) stop("\nxreg must have the same number of columns as the mean regressors in the model.")
        if (any(!is.finite(xreg))) stop("\nNA/NaN/Inf values found in xreg.")
        xreg <- coredata(xreg)
    } else {
        xreg <- matrix(0, nrow = h, ncol = ncols)
    }
    return(xreg)
}

.validate_sim_initv <- function(initv)
{
    if (any(!is.finite(initv)) || any(initv <= 0)) {
        stop("\nthe implied initial variance is not positive and finite; this can happen when the variance intercept is fixed at zero (e.g. ewma), which has no positive implied starting variance - supply a positive var_init.")
    }
    return(invisible(initv))
}

.simulate_garch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, xreg = NULL, burn = 0, ...)
{
    if (!is.null(seed)) set.seed(seed)
    extra_args <- list(...)
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    # combined pre-sample/burn-in length spanning both the GARCH order and
    # the ARMA order (mirroring rugarch's combined maxOrder, used uniformly
    # for both the variance and mean simulation recursions - see
    # .garchsimvec()/.armasimvec() in src/simulation.cpp).
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h + burn.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x (h + burn).")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }

    if (is.null(var_init)) {
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        p <- sum(alpha) + sum(beta)
        initv <- numerator/(1 - p)
    } else {
        initv <- var_init
    }
    .validate_sim_initv(initv)

    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }

    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_sqr_sim[,seq_len(maxpq)] <- matrix(initv, ncol = maxpq, nrow = nsim, byrow = TRUE)
        sigma_sim[,seq_len(maxpq)] <- matrix(sqrt(initv), ncol = maxpq, nrow = nsim, byrow = TRUE)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- as.integer(object$model$order)

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
    } else {
        # init column k feeds ARCH lag k, whose first lookback is pre-sample
        # column maxpq - k + 1 (column maxpq is the most recent pre-sample
        # period), so the chronological block must be reversed
        init <- epsilon[, rev(seq_len(maxpq)), drop = FALSE]^2
    }
    simc <- .garchsimvec(z = z, epsilon = epsilon, sigma_sqr_sim = sigma_sqr_sim, variance_intercept = variance_intercept, order = order, init = init, alpha = alpha, beta = beta, mu = mu, presample = maxpq)
    sigma <- simc$sigma
    series <- .arma_simulate_overlay(object, simc, mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_egarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, xreg = NULL, burn = 0, ...)
{
    if (!is.null(seed)) set.seed(seed)
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    extra_args <- list(...)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    kappa <- egarch_moment(distribution = distribution, skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta)
        initv <- (omega + mean_vreg)/(1 - p)
    } else {
        initv <- log(var_init)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        # innov_init is documented as applying identically to every sample
        # path, so it must be broadcast row-wise rather than assigned column-major
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }

    sigma_sim <- sigma_log_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_log_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(exp(initv))
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }

    order <- as.integer(object$model$order)

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
    } else {
        init <- (abs(z[, rev(seq_len(maxpq)), drop = FALSE]) - kappa)
    }

    simc <- .egarchsimvec(z = z, sigma_log_sim = sigma_log_sim, variance_intercept = variance_intercept, init = init, alpha = alpha, gamma = gamma, beta = beta, kappa = kappa, mu = mu, order = order, presample = maxpq)
    sigma <- simc$sigma
    series <- .arma_simulate_overlay(object, simc, mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_aparch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, xreg = NULL, burn = 0, ...)
{
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    if (!is.null(seed)) set.seed(seed)
    extra_args <- list(...)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }

    kappa <- aparch_moment_v(distribution = distribution, gamma = gamma, delta = delta,
                             skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta) + sum(alpha * kappa)
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - p)
    } else {
        initv <- var_init^(delta/2)
    }
    .validate_sim_initv(initv)
    order <- as.integer(object$model$order)
    sigma_sim <- sigma_power_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_power_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- initv^(1/delta)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]

        if (!is.null(extra_args$arch_initial)) {
            init <- extra_args$arch_initial
            init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
        } else {
            if (is.null(innov_init)) {
                init <- kappa * (initv^(delta/2))
                init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
            } else {
                # init column k feeds ARCH lag k: pre-sample column maxpq - k + 1
                # with gamma_k. Columns past order[1] are never read by the
                # recursion; the index is clamped to [1, order[1]] so that those
                # unread columns still hold a defined value rather than NA (an
                # order[1] of zero would otherwise index gamma with 0).
                e_pre <- epsilon[, rev(seq_len(maxpq)), drop = FALSE]
                g_idx <- pmax(pmin(seq_len(maxpq), order[1]), 1L)
                g_lag <- matrix(gamma[g_idx], ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
                init <- (abs(e_pre) - g_lag * e_pre)^delta
            }
        }
    }

    simc <- .aparchsimvec(epsilon = epsilon, sigma_power_sim = sigma_power_sim, z = z, variance_intercept = variance_intercept,
                          init = init, alpha = alpha, gamma = gamma, beta = beta, delta = delta, mu = mu, order = order, presample = maxpq)
    sigma <- simc$sigma
    series <- .arma_simulate_overlay(object, simc, mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_gjrgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, xreg = NULL, burn = 0, ...)
{
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    if (!is.null(seed)) set.seed(seed)
    extra_args <- list(...)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }
    kappa <- gjrgarch_moment(distribution = distribution, skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta) + sum(alpha) + sum(gamma * kappa)
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - p)
    } else {
        initv <- var_init
    }
    .validate_sim_initv(initv)

    sigma_sim <- sigma_squared_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_squared_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(initv)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
    } else {
        init <- (epsilon[, rev(seq_len(maxpq)), drop = FALSE]^2 * kappa)
    }

    order <- as.integer(object$model$order)
    simc <- .gjrsimvec(epsilon = epsilon, sigma_sqr_sim = sigma_squared_sim, z = z, variance_intercept = variance_intercept, order = order, init = init, alpha = alpha, gamma = gamma, beta = beta, mu = mu, presample = maxpq)

    sigma <- simc$sigma
    series <- .arma_simulate_overlay(object, simc, mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)

    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series, epsilon = epsilon)
    return(out)
}

.simulate_fgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, xreg = NULL, burn = 0, ...)
{
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    if (!is.null(seed)) set.seed(seed)
    extra_args <- list(...)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    eta <- object$parmatrix[group == "eta"]$value
    beta <- object$parmatrix[group == "beta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }

    kappa <- fgarch_moment_v(distribution = distribution, gamma = gamma, eta = eta, delta = delta,
                             skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta) + sum(alpha * kappa)
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - p)
    } else {
        initv <- var_init^(delta/2)
    }
    .validate_sim_initv(initv)

    sigma_sim <- sigma_power_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_power_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- initv^(1/delta)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }

    if (maxpq > 0) {
        if (!is.null(extra_args$arch_initial)) {
            init <- extra_args$arch_initial
            init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
        } else {
            if (is.null(innov_init)) {
                init <- kappa
                init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
            } else {
                # see .simulate_aparch() for the reversal and the index clamp;
                # note order is not yet in scope at this point in this function
                z_pre <- z[, rev(seq_len(maxpq)), drop = FALSE]
                g_idx <- pmax(pmin(seq_len(maxpq), object$model$order[1]), 1L)
                e_lag <- matrix(eta[g_idx], ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
                g_lag <- matrix(gamma[g_idx], ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
                d_pre <- z_pre - e_lag
                init <- (abs(d_pre) - g_lag * d_pre)^delta
            }
        }
    }

    order <- as.integer(object$model$order)

    simc <- .fgarchsimvec(epsilon = epsilon, sigma_power_sim = sigma_power_sim, z = z, variance_intercept = variance_intercept, init = init,
                          alpha = alpha, gamma = gamma, eta = eta, beta = beta, delta = delta, mu = mu, order = order, presample = maxpq)
    sigma <- simc$sigma
    series <- .arma_simulate_overlay(object, simc, mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_cgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, xreg = NULL, burn = 0, ...)
{
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    if (!is.null(seed)) set.seed(seed)
    extra_args <- list(...)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    rho <- object$parmatrix[group == "rho"]$value
    phi <- object$parmatrix[group == "phi"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }

    if (is.null(var_init)) {
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - rho)
        initq <- numerator/(1 - rho)
        if (object$vreg$multiplicative) {
            initv <- exp(initv)
            initq <- exp(initq)
        }
    } else {
        if (!is.matrix(var_init)) {
            if (length(var_init) != maxpq) stop(paste0("\nvar_init must be of length ", maxpq))
            initv <- var_init
            initq <- var_init
        } else {
            if (nrow(var_init) != maxpq) stop(paste0("\nvar_init must have ", maxpq, " rows"))
            initq <- var_init[,1]
            initv <- var_init[,2]
        }
    }
    .validate_sim_initv(initv)
    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    permanent_component_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    transitory_component_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        permanent_component_sim[,seq_len(maxpq)] <- initq
        sigma_sqr_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(initv)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- as.integer(object$model$order)
    for (i in (maxpq + 1):(h + maxpq)) {
        permanent_component_sim[,i] <- variance_intercept[i] + rho * permanent_component_sim[,i - 1] + phi * (epsilon[,i - 1]^2 - sigma_sqr_sim[,i - 1])
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                transitory_component_sim[,i] <- transitory_component_sim[,i] + alpha[j] * (epsilon[,i - j]^2 - sigma_sqr_sim[,i - j]) + alpha[j] * transitory_component_sim[,i - j]
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                transitory_component_sim[,i] <- transitory_component_sim[,i] + beta[j] * transitory_component_sim[,i - j]
            }
        }
        sigma_sqr_sim[,i] = transitory_component_sim[,i] + permanent_component_sim[,i]
        sigma_sim[,i] <- sqrt(sigma_sqr_sim[,i])
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    # return
    sigma <- sigma_sim
    permanent_component <- permanent_component_sim
    transitory_component <- transitory_component_sim
    series <- .arma_simulate_overlay(object, list(series = series_sim, epsilon = epsilon), mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        permanent_component <- permanent_component[,-seq_len(maxpq + burn), drop = FALSE]
        transitory_component <- transitory_component[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    class(permanent_component) <- "tsmodel.distribution"
    class(transitory_component) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    attr(permanent_component, "date_class") <- "numeric"
    attr(transitory_component, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series, transitory_component = transitory_component, permanent_component = permanent_component)
    return(out)
}


.simulate_igarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, xreg = NULL, burn = 0, ...)
{
    h <- h + burn
    xreg <- .process_simulation_xreg(object, xreg, h)
    if (!is.null(seed)) set.seed(seed)
    extra_args <- list(...)
    parameter <- group <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    maxpq <- max(object$model$order, arma_order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h + burn.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x (h + burn).")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }

    if (is.null(var_init)) {
        # set initiale variance close to integrated
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - 0.999)
    } else {
        initv <- var_init
    }
    .validate_sim_initv(initv)

    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(garch order, arma order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }

    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_sqr_sim[,seq_len(maxpq)] <- matrix(initv, ncol = maxpq, nrow = nsim, byrow = TRUE)
        sigma_sim[,seq_len(maxpq)] <- matrix(sqrt(initv), ncol = maxpq, nrow = nsim, byrow = TRUE)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- as.integer(object$model$order)

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        init <- .expand_arch_initial(init, maxpq, nrow(epsilon))
    } else {
        init <- epsilon[, rev(seq_len(maxpq)), drop = FALSE]^2
    }
    simc <- .garchsimvec(z = z, epsilon = epsilon, sigma_sqr_sim = sigma_sqr_sim, variance_intercept = variance_intercept, order = order, init = init, alpha = alpha, beta = beta, mu = mu, presample = maxpq)
    sigma <- simc$sigma
    series <- .arma_simulate_overlay(object, simc, mu, maxpq, arma_order, nsim, extra_args, xreg = xreg)
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}
