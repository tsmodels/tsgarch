# Stage B: regressor-aware filtering / prediction / simulation / backtesting.
# helper-global.R provides `y` (xts of dmbp). These tests use a 2-regressor
# arma(1,1)-garch(1,1) fit under both xreg_type conventions.

.make_armax_data <- function(n = 1200, seed = 42, m = 2) {
    set.seed(seed)
    dates <- as.Date(seq_len(n), origin = "1970-01-01")
    X <- xts(matrix(rnorm(n * m), ncol = m), dates)
    colnames(X) <- paste0("x", 1:m)
    yy <- xts(as.numeric(y[1:n, 1]) + as.numeric(coredata(X) %*% c(1.5, -0.8)[1:m]), dates)
    list(y = yy, x = X, dates = dates)
}

.fit_armax <- function(dat, xreg_type, n = NROW(dat$y), arma = c(1,1)) {
    spec <- garch_modelspec(dat$y[1:n], constant = TRUE, model = "garch", order = c(1,1),
                            arma = arma, xreg = dat$x[1:n], xreg_type = xreg_type)
    estimate(spec)
}

# hand-rolled point-forecast recursion for the conditional mean, mirroring
# the two conventions exactly
.hand_mean_forecast <- function(mod, newx, h) {
    group <- NULL
    ac <- arma_coefficients(mod)
    ar <- as.numeric(ac$ar)
    ma <- as.numeric(ac$ma)
    ar_order <- length(ar)
    ma_order <- length(ma)
    maxpq <- max(ar_order, ma_order)
    mu <- mod$parmatrix[group == "mu"]$value
    tau <- mod$parmatrix[group == "tau"]$value
    armax <- identical(mod$spec$xreg$xreg_type, "armax")
    xtau_f <- as.numeric(as.matrix(newx) %*% tau)
    if (maxpq == 0) return(mu + xtau_f)
    y_hist <- tail(as.numeric(mod$spec$target$y), maxpq)
    eps_hist <- tail(as.numeric(residuals(mod)), maxpq)
    xtau_hist <- tail(as.numeric(mod$spec$xreg$xreg %*% tau), maxpq)
    x <- if (armax) c(y_hist - mu, rep(0, h)) else c(y_hist - mu - xtau_hist, rep(0, h))
    eps <- c(eps_hist, rep(0, h))
    for (i in 1:h) {
        idx <- maxpq + i
        val <- if (armax) xtau_f[i] else 0
        if (ar_order > 0) for (j in 1:ar_order) val <- val + ar[j] * x[idx - j]
        if (ma_order > 0) for (j in 1:ma_order) if (i - j <= 0) val <- val + ma[j] * eps[idx - j]
        x[idx] <- val
    }
    if (armax) mu + x[(maxpq + 1):(maxpq + h)] else mu + xtau_f + x[(maxpq + 1):(maxpq + h)]
}

test_that("tsfilter: incremental and full-sample filtering agree (both conventions)", {
    dat <- .make_armax_data()
    N <- 1000
    for (xt in c("arma_errors","armax")) {
        mod <- .fit_armax(dat, xt, n = N)
        y_new <- dat$y[(N + 1):1200]
        x_new <- dat$x[(N + 1):1200]
        f_all <- tsfilter(mod, y = y_new, newxreg = x_new)
        f_inc <- mod
        for (k in seq_len(NROW(y_new))) {
            f_inc <- tsfilter(f_inc, y = y_new[k], newxreg = x_new[k])
        }
        # full-sample reference: spec with the fitted parameters held fixed,
        # filtered over all 1200 observations in one pass
        spec_full <- garch_modelspec(dat$y, constant = TRUE, model = "garch", order = c(1,1),
                                     arma = c(1,1), xreg = dat$x, xreg_type = xt)
        spec_full$parmatrix$value <- f_all$parmatrix$value
        f_full <- tsfilter(spec_full)
        n_new <- NROW(y_new)
        expect_equal(as.numeric(fitted(f_inc)), as.numeric(fitted(f_full)), tolerance = 1e-8,
                     info = paste("fitted", xt))
        expect_equal(as.numeric(residuals(f_inc)), as.numeric(residuals(f_full)), tolerance = 1e-8,
                     info = paste("residuals", xt))
        # incremental (one-at-a-time) vs batch filtering take the same code
        # path and must agree over the whole sample
        expect_equal(as.numeric(sigma(f_inc)), as.numeric(sigma(f_all)), tolerance = 1e-8,
                     info = paste("sigma inc-vs-batch", xt))
        expect_equal(as.numeric(fitted(f_all)), as.numeric(fitted(f_full)), tolerance = 1e-8,
                     info = paste("fitted batch", xt))
        # sigma over the newly filtered segment: agreement with the
        # full-sample spec filter is exact here. The early in-sample region
        # differs, but not because the two code paths disagree - on an
        # identical sample they are bit-identical. The cause is that the
        # recursion seed (init = "unconditional", i.e. the mean of the
        # squared residuals) is computed over whatever sample is being
        # filtered, so seeding from 1000 observations and appending 200 is
        # not the same as filtering 1200 in one pass. The resulting
        # difference in sigma^2 decays at exactly beta per period (the ARCH
        # term is unaffected, since the residuals are identical), so it is a
        # geometric transient that is numerically dead long before the
        # filtered segment. Unrelated to xreg, and present for xreg-free
        # models too.
        expect_equal(tail(as.numeric(sigma(f_inc)), n_new), tail(as.numeric(sigma(f_full)), n_new),
                     tolerance = 1e-8, info = paste("sigma new segment", xt))
        expect_equal(tail(as.numeric(sigma(f_all)), n_new), tail(as.numeric(sigma(f_full)), n_new),
                     tolerance = 1e-8, info = paste("sigma new segment batch", xt))
    }
})

test_that("predict: point forecasts match hand-rolled recursion (both conventions, h = 1 and 5)", {
    dat <- .make_armax_data()
    for (xt in c("arma_errors","armax")) {
        mod <- .fit_armax(dat, xt)
        for (h in c(1, 5)) {
            newx <- xts(matrix(rnorm(h * 2), ncol = 2), .forecast_dates(NULL, h, "days", tail(dat$dates, 1)))
            colnames(newx) <- colnames(dat$x)
            p <- predict(mod, h = h, newxreg = newx, nsim = 0)
            expect_equal(as.numeric(p$mean), .hand_mean_forecast(mod, newx, h),
                         tolerance = 1e-8, info = paste(xt, "h =", h))
        }
    }
})

test_that("predict: arma = c(0,0) with xreg gives mean = mu + newxreg %*% tau", {
    group <- NULL
    dat <- .make_armax_data()
    for (xt in c("arma_errors","armax")) {
        mod <- .fit_armax(dat, xt, arma = c(0,0))
        h <- 4
        newx <- matrix(rnorm(h * 2), ncol = 2)
        p <- predict(mod, h = h, newxreg = newx, nsim = 0)
        mu <- mod$parmatrix[group == "mu"]$value
        tau <- mod$parmatrix[group == "tau"]$value
        expect_equal(as.numeric(p$mean), as.numeric(mu + newx %*% tau), tolerance = 1e-10,
                     info = xt)
    }
})

# hand-rolled deterministic simulation matching .garchsimvec + .armasimvec;
# igarch = TRUE reproduces .simulate_igarch()'s initv = numerator/(1 - 0.999)
.hand_simulate <- function(spec, innov, xreg = NULL, igarch = FALSE) {
    group <- parameter <- NULL
    arma_order <- spec$model$arma
    maxpq <- max(spec$model$order, arma_order)
    h <- length(innov)
    mu <- spec$parmatrix[parameter == "mu"]$value
    omega <- spec$parmatrix[parameter == "omega"]$value
    alpha <- spec$parmatrix[group == "alpha"]$value
    beta <- spec$parmatrix[group == "beta"]$value
    tau <- spec$parmatrix[group == "tau"]$value
    ac <- arma_coefficients(spec)
    ar <- as.numeric(ac$ar)
    ma <- as.numeric(ac$ma)
    armax <- identical(spec$xreg$xreg_type, "armax")
    T <- maxpq + h
    z <- c(rep(1, maxpq), innov)
    # variance intercept mirrors .simulate_garch()/.simulate_igarch(): omega
    # (or exp(omega) under a multiplicative specification with vreg)
    vint <- if (spec$vreg$multiplicative) exp(omega) else omega
    initv <- if (igarch) vint/(1 - 0.999) else vint/(1 - sum(alpha) - sum(beta))
    sig2 <- c(rep(initv, maxpq), rep(0, h))
    eps <- numeric(T)
    eps[seq_len(maxpq)] <- z[seq_len(maxpq)] * sqrt(initv)
    for (i in (maxpq + 1):T) {
        sig2[i] <- vint + alpha[1] * eps[i - 1]^2 + beta[1] * sig2[i - 1]
        eps[i] <- z[i] * sqrt(sig2[i])
    }
    if (is.null(xreg)) xreg <- matrix(0, nrow = h, ncol = NCOL(spec$xreg$xreg))
    xtau <- c(as.numeric(tail(spec$xreg$xreg %*% tau, maxpq)),
              as.numeric(as.matrix(xreg) %*% tau))
    series <- numeric(T)
    series[seq_len(maxpq)] <- if (armax) mu else mu + xtau[seq_len(maxpq)]
    ma_eps <- eps
    ma_eps[seq_len(maxpq)] <- 0
    for (i in (maxpq + 1):T) {
        v <- mu + xtau[i]
        if (length(ar) > 0) {
            for (j in seq_along(ar)) {
                dev <- series[i - j] - mu - if (armax) 0 else xtau[i - j]
                v <- v + ar[j] * dev
            }
        }
        if (length(ma) > 0) for (j in seq_along(ma)) v <- v + ma[j] * ma_eps[i - j]
        series[i] <- v + eps[i]
    }
    series[(maxpq + 1):T]
}

test_that("simulate: fixed innovations reproduce hand-rolled recursion (both conventions)", {
    dat <- .make_armax_data()
    for (xt in c("arma_errors","armax")) {
        mod <- .fit_armax(dat, xt)
        # the spec embedded in the estimate has parmatrix = NULL by design;
        # rebuild a standalone spec carrying the fitted parameter values
        spec <- garch_modelspec(dat$y, constant = TRUE, model = "garch", order = c(1,1),
                                arma = c(1,1), xreg = dat$x, xreg_type = xt)
        spec$parmatrix <- data.table::copy(mod$parmatrix)
        h <- 6
        set.seed(7)
        innov <- matrix(rnorm(h), nrow = 1)
        xnew <- matrix(rnorm(h * 2), ncol = 2)
        s <- simulate(spec, nsim = 1, h = h, innov = innov, xreg = xnew)
        expect_equal(as.numeric(s$series[1,]), .hand_simulate(spec, as.numeric(innov), xnew),
                     tolerance = 1e-10, info = xt)
    }
})

test_that("predict: parametric distribution is centered on the analytic mean", {
    dat <- .make_armax_data()
    for (xt in c("arma_errors","armax")) {
        mod <- .fit_armax(dat, xt)
        h <- 5
        newx <- xts(matrix(rnorm(h * 2), ncol = 2), .forecast_dates(NULL, h, "days", tail(dat$dates, 1)))
        p <- predict(mod, h = h, newxreg = newx, nsim = 4000, sim_method = "parametric", seed = 11)
        center <- colMeans(as.matrix(p$distribution))
        expect_equal(as.numeric(center), as.numeric(p$mean), tolerance = 0.02,
                     info = xt)
    }
})

test_that("validation: newxreg/xreg misuse errors and NULL warns to zero", {
    dat <- .make_armax_data()
    N <- 1000
    mod <- .fit_armax(dat, "armax", n = N)
    # a spec carrying xreg (simulate() only needs the spec + parmatrix; the
    # spec embedded in an estimate object has parmatrix = NULL by design)
    specx <- garch_modelspec(dat$y[1:N], constant = TRUE, model = "garch", order = c(1,1),
                             arma = c(1,1), xreg = dat$x[1:N], xreg_type = "armax")
    specx$parmatrix <- data.table::copy(mod$parmatrix)
    h <- 3
    newx <- matrix(rnorm(h * 2), ncol = 2)
    # NULL warns and is equivalent to zeros
    expect_warning(p_null <- predict(mod, h = h, newxreg = NULL, nsim = 0), "newxreg")
    p_zero <- predict(mod, h = h, newxreg = matrix(0, h, 2), nsim = 0)
    expect_equal(as.numeric(p_null$mean), as.numeric(p_zero$mean))
    expect_warning(simulate(specx, h = h, xreg = NULL), "xreg")
    # with fixed innovations the NULL-vs-zero paths are identical
    innov <- matrix(rnorm(h), nrow = 1)
    s_null2 <- suppressWarnings(simulate(specx, h = h, innov = innov, xreg = NULL))
    s_zero2 <- simulate(specx, h = h, innov = innov, xreg = matrix(0, h, 2))
    expect_equal(as.numeric(s_null2$series), as.numeric(s_zero2$series))
    y_new <- dat$y[(N + 1):(N + h)]
    expect_warning(f_null <- tsfilter(mod, y = y_new, newxreg = NULL), "newxreg")
    f_zero <- tsfilter(.fit_armax(dat, "armax", n = N), y = y_new, newxreg = matrix(0, h, 2))
    expect_equal(as.numeric(sigma(f_null)), as.numeric(sigma(f_zero)))
    # wrong ncol / nrow / non-finite
    expect_error(predict(mod, h = h, newxreg = matrix(0, h, 3), nsim = 0), "columns")
    expect_error(predict(mod, h = h, newxreg = matrix(0, h + 1, 2), nsim = 0), "rows")
    expect_error(predict(mod, h = h, newxreg = matrix(NA_real_, h, 2), nsim = 0), "NA/NaN/Inf")
    expect_error(predict(mod, h = h, newxreg = matrix(Inf, h, 2), nsim = 0), "NA/NaN/Inf")
    expect_error(simulate(specx, h = h, xreg = matrix(0, h, 3)), "columns")
    expect_error(simulate(specx, h = h, xreg = matrix(0, h + 1, 2)), "rows")
    expect_error(simulate(specx, h = h, xreg = matrix(NaN, h, 2)), "NA/NaN/Inf")
    expect_error(tsfilter(.fit_armax(dat, "armax", n = N), y = y_new, newxreg = matrix(0, h, 3)), "columns")
    expect_error(tsfilter(.fit_armax(dat, "armax", n = N), y = y_new, newxreg = matrix(0, h + 1, 2)), "rows")
    expect_error(tsfilter(.fit_armax(dat, "armax", n = N), y = y_new, newxreg = matrix(-Inf, h, 2)), "NA/NaN/Inf")
})

test_that(".spec2newspec round-trips xreg and xreg_type; tsfilter on spec runs", {
    dat <- .make_armax_data()
    for (xt in c("arma_errors","armax")) {
        mod <- .fit_armax(dat, xt)
        spec2 <- tsgarch:::.spec2newspec(mod$spec)
        expect_true(spec2$xreg$include_xreg)
        expect_equal(spec2$xreg$xreg_type, xt)
        expect_equal(spec2$xreg$xreg, mod$spec$xreg$xreg)
        expect_equal(NROW(spec2$parmatrix[group == "tau"]), 2)
        spec2$parmatrix$value <- mod$parmatrix$value
        f <- tsfilter(spec2)
        expect_s3_class(f, "tsgarch.estimate")
        expect_equal(as.numeric(fitted(f)), as.numeric(fitted(mod)), tolerance = 1e-8, info = xt)
    }
    # NULL-safe fallback for objects without an xreg slot
    spec_old <- mod$spec
    spec_old$xreg <- NULL
    spec3 <- tsgarch:::.spec2newspec(spec_old)
    expect_false(spec3$xreg$include_xreg)
    expect_equal(spec3$xreg$xreg_type, "arma_errors")
})

test_that("simulate: igarch/ewma honor the ARMA overlay (previously ignored)", {
    dat <- .make_armax_data()
    h <- 6
    set.seed(7)
    innov <- matrix(rnorm(h), nrow = 1)
    for (m in c("igarch","ewma")) {
        spec <- garch_modelspec(dat$y, constant = TRUE, model = m, order = c(1,1), arma = c(1,1))
        # ewma fixes omega at 0 (variance intercept collapses to zero); give
        # it a small positive value so the simulated sigma path is exercised
        spec$parmatrix[parameter == "omega", value := 0.02]
        s <- simulate(spec, nsim = 1, h = h, innov = innov)
        expect_equal(as.numeric(s$series[1,]), .hand_simulate(spec, as.numeric(innov), igarch = TRUE),
                     tolerance = 1e-10, info = m)
        # the overlay is genuinely active: the series is not just mu + eps
        mu <- spec$parmatrix[parameter == "mu"]$value
        eps_path <- mu + as.numeric(innov) * as.numeric(s$sigma[1,])
        expect_gt(max(abs(as.numeric(s$series[1,]) - eps_path)), 1e-4)
    }
})

test_that("simulate: igarch with xreg honors both conventions", {
    dat <- .make_armax_data()
    h <- 6
    set.seed(7)
    innov <- matrix(rnorm(h), nrow = 1)
    xnew <- matrix(rnorm(h * 2), ncol = 2)
    for (xt in c("arma_errors","armax")) {
        spec <- garch_modelspec(dat$y, constant = TRUE, model = "igarch", order = c(1,1),
                                arma = c(1,1), xreg = dat$x, xreg_type = xt)
        s <- simulate(spec, nsim = 1, h = h, innov = innov, xreg = xnew)
        expect_equal(as.numeric(s$series[1,]),
                     .hand_simulate(spec, as.numeric(innov), xnew, igarch = TRUE),
                     tolerance = 1e-10, info = xt)
    }
})

test_that("simulate: igarch arma = c(0,0) with no xreg is unchanged (mu + eps)", {
    dat <- .make_armax_data()
    h <- 6
    set.seed(7)
    innov <- matrix(rnorm(h), nrow = 1)
    spec <- garch_modelspec(dat$y, constant = TRUE, model = "igarch", order = c(1,1), arma = c(0,0))
    s <- simulate(spec, nsim = 1, h = h, innov = innov)
    mu <- spec$parmatrix[parameter == "mu"]$value
    # the old path returned series = mu + epsilon exactly; pin that
    expect_equal(as.numeric(s$series[1,]), mu + as.numeric(innov) * as.numeric(s$sigma[1,]),
                 tolerance = 1e-12)
})

test_that("simulate: spec without parmatrix values errors cleanly (no abort)", {
    mod <- .fit_armax(.make_armax_data(n = 400), "armax")
    # the spec embedded in an estimate carries parmatrix = NULL; assigning
    # NULL to the value column used to drop it and abort the session
    expect_error(simulate(mod$spec, h = 2), "parmatrix")
})

test_that("backtest: first-window mu matches manual refit + predict (both conventions)", {
    dat <- .make_armax_data(n = 400)
    for (xt in c("arma_errors","armax")) {
        spec <- garch_modelspec(dat$y, constant = TRUE, model = "garch", order = c(1,1),
                                arma = c(1,1), xreg = dat$x, xreg_type = xt)
        start <- 360
        bt <- tsbacktest(spec, start = start, end = 375, h = 1, estimate_every = 1, rolling = FALSE)
        first_date <- bt$table$estimation_date[1]
        fc_date <- bt$table$forecast_date[1]
        y_train <- dat$y[paste0("/", first_date)]
        x_train <- dat$x[index(y_train)]
        mod <- estimate(garch_modelspec(y_train, constant = TRUE, model = "garch", order = c(1,1),
                                        arma = c(1,1), xreg = x_train, xreg_type = xt))
        p <- predict(mod, h = 1, newxreg = dat$x[fc_date], forc_dates = fc_date, nsim = 0)
        expect_equal(bt$table$mu[1], as.numeric(p$mean), tolerance = 1e-8, info = xt)
    }
})
