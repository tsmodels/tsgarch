test_that("armax: xreg = NULL changes nothing (spec shape and fit unchanged)", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1))
    # tau dummy row is always present (keeps parmatrix aligned with the
    # always >= 1 column x matrix passed to TMB) but never estimated
    expect_equal(NROW(spec$parmatrix[group == "tau"]), 1)
    expect_equal(spec$parmatrix[group == "tau"]$estimate, 0)
    expect_equal(spec$parmatrix[group == "tau"]$parameter, "tau1")
    expect_false(spec$xreg$include_xreg)
    expect_equal(spec$xreg$xreg_type, "arma_errors")
    expect_length(spec$model_options, 9)
    # a spec with explicit xreg = NULL is identical to one where it is omitted
    spec2 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1),
                             xreg = NULL, xreg_type = "arma_errors")
    mod <- estimate(spec); mod2 <- estimate(spec2)
    expect_equal(as.numeric(logLik(mod)), as.numeric(logLik(mod2)))
    # the no-xreg logLik values are pinned by test-estimation.R / test-arma.R,
    # which continue to pass unchanged (regression guard for this feature)
    speca <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    moda <- estimate(speca)
    expect_true(is.finite(as.numeric(logLik(moda))))
})

# helper: simulate regression + garch(1,1) errors
.sim_armax_errors <- function(mu = 0.5, tau = 2, alpha = 0.1, beta = 0.85, omega = 0.05,
                              phi = 0, n = 2500, seed = 1, m = 1) {
    set.seed(seed)
    X <- matrix(rnorm(n * m), ncol = m)
    if (m == 1) tau_vec <- tau else tau_vec <- tau
    sigma2 <- eps <- w <- yy <- numeric(n)
    sigma2[1] <- omega/(1 - alpha - beta)
    for (i in 2:n) {
        sigma2[i] <- omega + alpha * eps[i - 1]^2 + beta * sigma2[i - 1]
        eps[i] <- rnorm(1) * sqrt(sigma2[i])
        w[i] <- phi * w[i - 1] + eps[i]
        yy[i] <- mu + sum(X[i,] * tau_vec) + w[i]
    }
    list(y = xts(yy, as.Date(seq_along(yy), origin = "1970-01-01")),
         x = xts(X, as.Date(seq_along(yy), origin = "1970-01-01")))
}

test_that("armax: pure regression + GARCH(1,1) recovers mu and tau", {
    d <- .sim_armax_errors(mu = 0.5, tau = 2, n = 2500, seed = 5)
    spec <- garch_modelspec(d$y, model = "garch", constant = TRUE, order = c(1,1),
                            arma = c(0,0), xreg = d$x)
    mod <- estimate(spec)
    expect_equal(mod$parmatrix[parameter == "mu"]$value, 0.5, tolerance = 0.05)
    expect_equal(mod$parmatrix[parameter == "tau1"]$value, 2, tolerance = 0.05)
})

test_that("armax: arma_errors and armax coincide iff ar_order == 0", {
    d <- .sim_armax_errors(mu = 0.5, tau = 2, n = 2500, seed = 6)
    spec_ae <- garch_modelspec(d$y, model = "garch", constant = TRUE, order = c(1,1),
                               arma = c(0,1), xreg = d$x, xreg_type = "arma_errors")
    spec_ax <- garch_modelspec(d$y, model = "garch", constant = TRUE, order = c(1,1),
                               arma = c(0,1), xreg = d$x, xreg_type = "armax")
    mod_ae <- estimate(spec_ae); mod_ax <- estimate(spec_ax)
    expect_equal(as.numeric(logLik(mod_ae)), as.numeric(logLik(mod_ax)), tolerance = 1e-6)
    expect_equal(mod_ae$parmatrix[parameter == "tau1"]$value,
                 mod_ax$parmatrix[parameter == "tau1"]$value, tolerance = 1e-6)
    # with a nonzero AR order the two conventions are different models
    # (data must have genuine AR dynamics for the difference to show, so
    # phi = 0.6 here rather than the iid-errors default)
    d2 <- .sim_armax_errors(mu = 0.5, tau = 2, phi = 0.6, n = 2500, seed = 61)
    spec_ae1 <- garch_modelspec(d2$y, model = "garch", constant = TRUE, order = c(1,1),
                                arma = c(1,0), xreg = d2$x, xreg_type = "arma_errors")
    spec_ax1 <- garch_modelspec(d2$y, model = "garch", constant = TRUE, order = c(1,1),
                                arma = c(1,0), xreg = d2$x, xreg_type = "armax")
    mod_ae1 <- estimate(spec_ae1); mod_ax1 <- estimate(spec_ax1)
    expect_false(isTRUE(all.equal(as.numeric(logLik(mod_ae1)), as.numeric(logLik(mod_ax1)), tolerance = 1e-4)))
    expect_false(isTRUE(all.equal(mod_ae1$parmatrix[parameter == "tau1"]$value,
                                  mod_ax1$parmatrix[parameter == "tau1"]$value, tolerance = 1e-4)))
})

test_that("armax: arma_errors matches stats::arima(xreg=) under constant variance", {
    set.seed(7)
    n <- 2500
    X <- rnorm(n)
    e <- arima.sim(n = n, list(ar = 0.6))
    yy <- 0.5 + 2 * X + as.numeric(e)
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    xs <- xts(matrix(X, ncol = 1), index(ys))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0),
                            arma = c(1,0), xreg = xs)
    mod <- estimate(spec)
    fit <- arima(yy, order = c(1,0,0), xreg = X)
    ac <- arma_coefficients(mod)
    expect_equal(mod$parmatrix[parameter == "mu"]$value, unname(coef(fit)["intercept"]), tolerance = 5e-3)
    expect_equal(unname(ac$ar), unname(coef(fit)["ar1"]), tolerance = 1e-3)
    expect_equal(mod$parmatrix[parameter == "tau1"]$value, unname(coef(fit)["X"]), tolerance = 1e-3)
})

test_that("armax: xreg_type = 'armax' recovers its own DGP (impact effect)", {
    # persistent regressor so the lagged-regressor term (by which the two
    # conventions differ) is genuinely informative
    set.seed(8)
    n <- 3000
    X <- as.numeric(arima.sim(n = n, list(ar = 0.9)))
    mu_true <- 0.5; phi <- 0.6; tau_true <- 1.5
    eps <- rnorm(n)
    yy <- numeric(n)
    for (i in 2:n) yy[i] <- mu_true + phi * (yy[i - 1] - mu_true) + tau_true * X[i] + eps[i]
    yy <- yy[500:n]
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    xs <- xts(matrix(X[500:n], ncol = 1), index(ys))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0),
                            arma = c(1,0), xreg = xs, xreg_type = "armax")
    mod <- suppressWarnings(estimate(spec))
    ac <- arma_coefficients(mod)
    expect_equal(unname(ac$ar), phi, tolerance = 0.05)
    expect_equal(mod$parmatrix[parameter == "tau1"]$value, tau_true, tolerance = 0.05)
    # the conventions differ by a single -phi*tau*x_{t-1} term, so the
    # misspecified arma_errors fit is a different model: on this data it
    # absorbs the omitted lagged-regressor term mostly through a much larger
    # AR coefficient (ar: 0.61 vs 0.95, ll: -3610 vs -4636) rather than
    # through tau (1.49 vs 1.38 - the tau gap alone is small and
    # seed-dependent, so the AR coefficient and the likelihood are the
    # sharper discriminators)
    spec_ae <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0),
                               arma = c(1,0), xreg = xs, xreg_type = "arma_errors")
    mod_ae <- suppressWarnings(estimate(spec_ae))
    expect_true(abs(unname(arma_coefficients(mod_ae)$ar) - unname(ac$ar)) > 0.2)
    expect_true(abs(as.numeric(logLik(mod_ae)) - as.numeric(logLik(mod))) > 10)
    expect_false(isTRUE(all.equal(mod_ae$parmatrix[parameter == "tau1"]$value,
                                  mod$parmatrix[parameter == "tau1"]$value, tolerance = 1e-3)))
})

# hand-rolled conditional mean for both conventions
.handrolled_fitted <- function(mod, xmat) {
    ac <- arma_coefficients(mod)
    ar <- as.numeric(ac$ar); ma <- as.numeric(ac$ma)
    mu <- mod$parmatrix[parameter == "mu"]$value
    tau <- mod$parmatrix[group == "tau"]$value
    yv <- as.numeric(mod$spec$target$y_orig)
    xtau <- drop(xmat %*% tau)
    armax <- identical(mod$spec$xreg$xreg_type, "armax")
    n <- length(yv)
    # one zero pre-sample row, matching the cmodel(0) burn-in convention
    z <- c(0, if (armax) yv - mu else yv - mu - xtau)
    eps <- rep(0, n + 1)
    cond <- rep(mu, n + 1)
    for (t in 2:(n + 1)) {
        m_ <- 0
        if (length(ar) > 0) m_ <- m_ + ar[1] * z[t - 1]
        if (length(ma) > 0) m_ <- m_ + ma[1] * eps[t - 1]
        if (armax) m_ <- m_ + xtau[t - 1]
        eps[t] <- z[t] - m_
        cond[t] <- yv[t - 1] - eps[t]
    }
    cond[-1]
}

test_that("armax: fitted() equals the hand-rolled recursion and residuals() == y - fitted(), both conventions", {
    set.seed(9)
    n <- 1500
    X <- cbind(rnorm(n), rnorm(n))
    e <- as.numeric(arima.sim(n = n, list(ar = 0.5, ma = 0.3)))
    yy <- 0.3 + X %*% c(1.5, -0.7) + e
    ys <- xts(as.numeric(yy), as.Date(seq_along(yy), origin = "1970-01-01"))
    xs <- xts(X, index(ys))
    for (tp in c("arma_errors", "armax")) {
        spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(1,1),
                                arma = c(1,1), xreg = xs, xreg_type = tp)
        mod <- estimate(spec)
        expect_equal(as.numeric(residuals(mod)), as.numeric(mod$spec$target$y_orig) - as.numeric(fitted(mod)))
        expect_equal(as.numeric(fitted(mod)), .handrolled_fitted(mod, X), tolerance = 1e-8)
    }
})

test_that("armax: tau rows can be fixed (estimate = 0) and drop out of coef()", {
    d <- .sim_armax_errors(mu = 0.5, tau = c(1.5, -0.7), n = 2000, seed = 10, m = 2)
    spec <- garch_modelspec(d$y, model = "garch", constant = TRUE, order = c(1,1),
                            arma = c(0,0), xreg = d$x)
    spec$parmatrix[group == "tau", value := c(1.5, -0.7)]
    spec$parmatrix[group == "tau", estimate := 0]
    mod <- estimate(spec)
    expect_equal(mod$parmatrix[group == "tau"]$value, c(1.5, -0.7))
    expect_false(any(grepl("^tau", names(coef(mod)))))
    expect_true(is.finite(as.numeric(logLik(mod))))
})

test_that("armax: all 8 model flavours estimate with xreg under both conventions", {
    d <- .sim_armax_errors(mu = 0.5, tau = c(1.2, -0.5), n = 1500, seed = 11, m = 2)
    for (mdl in c("garch", "egarch", "gjrgarch", "aparch", "fgarch", "cgarch", "igarch", "ewma")) {
        for (tp in c("arma_errors", "armax")) {
            spec <- garch_modelspec(d$y, model = mdl, constant = TRUE, order = c(1,1),
                                    arma = c(1,1), xreg = d$x, xreg_type = tp)
            mod <- suppressWarnings(estimate(spec))
            expect_true(is.finite(as.numeric(logLik(mod))), info = paste(mdl, tp))
            expect_true(all(c("tau1","tau2") %in% names(coef(mod))), info = paste(mdl, tp))
        }
    }
})

test_that("armax: specification-time validation of xreg", {
    ys <- y[1:500,1]
    xs <- xts(matrix(rnorm(500), ncol = 1), index(ys))
    # rank-deficient design
    x_bad <- xts(cbind(xs, 2 * xs[,1]), index(ys))
    expect_error(garch_modelspec(ys, model = "garch", constant = FALSE, xreg = x_bad), "rank deficient")
    # constant column collinear with the estimated constant
    x_const <- xts(cbind(rep(1, 500), rnorm(500)), index(ys))
    expect_error(garch_modelspec(ys, model = "garch", constant = TRUE, xreg = x_const), "rank deficient")
    # wrong number of rows
    x_short <- xts(matrix(rnorm(400), ncol = 1), index(ys)[1:400])
    expect_error(garch_modelspec(ys, model = "garch", xreg = x_short), "same number of rows")
    # non-finite values
    x_inf <- xs; x_inf[10,1] <- Inf
    expect_error(garch_modelspec(ys, model = "garch", xreg = x_inf), "NA/NaN/Inf")
    x_na <- xs; x_na[10,1] <- NA
    expect_error(garch_modelspec(ys, model = "garch", xreg = x_na), "NA/NaN/Inf")
    # invalid xreg_type
    expect_error(garch_modelspec(ys, model = "garch", xreg = xs, xreg_type = "bogus"))
})
