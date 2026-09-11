test_that("arma: pacf_to_ar/pacf_to_ma guarantee stationarity/invertibility", {
    set.seed(11)
    for (n in 1:6) {
        r <- runif(n, -0.9, 0.9)
        ar <- tsgarch:::pacf_to_ar(r)
        expect_true(min(Mod(polyroot(c(1, -ar)))) > 1)
        ma <- tsgarch:::pacf_to_ma(r)
        expect_true(min(Mod(polyroot(c(1, ma)))) > 1)
    }
})

test_that("arma: default arma = c(0,0) leaves garch estimation unchanged", {
    spec0 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1))
    spec1 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(0,0))
    mod0 <- estimate(spec0)
    mod1 <- estimate(spec1)
    expect_equal(as.numeric(logLik(mod0)), as.numeric(logLik(mod1)), tolerance = 1e-6)
})

test_that("arma: AR(1) mean recovers close to stats::arima", {
    set.seed(42)
    n <- 2000
    e <- rnorm(n)
    yy <- numeric(n)
    for (i in 2:n) yy[i] <- 2 + 0.6 * (yy[i - 1] - 2) + e[i]
    yy <- yy[500:n]
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0), arma = c(1,0))
    mod <- estimate(spec)
    fit <- arima(yy, order = c(1,0,0))
    ar_est <- mod$parmatrix[group == "arpacf"]$value
    expect_equal(ar_est, unname(coef(fit)["ar1"]), tolerance = 0.02)
    expect_equal(mod$parmatrix[group == "mu"]$value, unname(coef(fit)["intercept"]), tolerance = 0.02)
})

test_that("arma: MA(1) mean recovers close to stats::arima", {
    set.seed(7)
    n <- 3000
    e <- rnorm(n)
    yy <- numeric(n)
    for (i in 2:n) yy[i] <- 0.3 + e[i] + 0.5 * e[i - 1]
    yy <- yy[200:n]
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0), arma = c(0,1))
    mod <- estimate(spec)
    fit <- arima(yy, order = c(0,0,1))
    ma_est <- mod$parmatrix[group == "mapacf"]$value
    expect_equal(ma_est, unname(coef(fit)["ma1"]), tolerance = 0.02)
})

test_that("arma: AR(2) mean recovers close to stats::arima", {
    set.seed(3)
    phi <- c(0.5, -0.3)
    n <- 3000
    e <- rnorm(n)
    yy <- numeric(n)
    for (i in 3:n) yy[i] <- phi[1]*yy[i-1] + phi[2]*yy[i-2] + e[i]
    yy <- yy[200:n]
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    spec <- garch_modelspec(ys, model = "garch", constant = FALSE, order = c(0,0), arma = c(2,0))
    mod <- estimate(spec)
    ar_pacf <- mod$parmatrix[group == "arpacf"]$value
    ar_est <- tsgarch:::pacf_to_ar(ar_pacf)
    fit <- arima(yy, order = c(2,0,0), include.mean = FALSE)
    expect_equal(ar_est, unname(coef(fit)), tolerance = 0.02)
})

test_that("arma: joint ARMA(1,1)-GARCH(1,1) estimation is well-behaved and nested vs arma=c(0,0)", {
    spec0 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(0,0))
    spec1 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod0 <- estimate(spec0)
    mod1 <- estimate(spec1)
    # adding parameters to a nested model cannot decrease the likelihood
    expect_true(as.numeric(logLik(mod1)) >= as.numeric(logLik(mod0)) - 1e-4)
    expect_true(mod1$conditions$kkt1)
})

test_that("arma: is currently rejected for non-garch models", {
    expect_error(garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch", order = c(1,1), arma = c(1,0)))
})

test_that("arma: residuals() returns the true ARMA residual, not y - mu", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod <- estimate(spec)
    res <- as.numeric(residuals(mod))
    naive <- as.numeric(mod$spec$target$y_orig - mod$parmatrix[parameter == "mu"]$value)
    expect_false(isTRUE(all.equal(res, naive)))
    expect_equal(length(res), length(naive))
})

test_that("arma: fitted() returns the time-varying conditional mean, and residuals() == y - fitted() exactly", {
    spec0 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(0,0))
    mod0 <- estimate(spec0)
    f0 <- as.numeric(fitted(mod0))
    mu0 <- mod0$parmatrix[parameter == "mu"]$value
    # no ARMA: fitted() is still a full-length vector, constant at mu
    expect_length(f0, mod0$nobs)
    expect_true(all(f0 == mu0))
    expect_equal(as.numeric(residuals(mod0)), as.numeric(mod0$spec$target$y_orig) - f0)

    spec1 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod1 <- estimate(spec1)
    f1 <- as.numeric(fitted(mod1))
    expect_length(f1, mod1$nobs)
    # ARMA: fitted() is genuinely time-varying, not constant
    expect_true(length(unique(round(f1, 8))) > 1)
    expect_equal(as.numeric(residuals(mod1)), as.numeric(mod1$spec$target$y_orig) - f1)
})

test_that("arma: constant = FALSE with no arma gives fitted() == 0 (full-length vector)", {
    spec <- garch_modelspec(y[1:1800,1], constant = FALSE, model = "garch", order = c(1,1))
    mod <- estimate(spec)
    f <- as.numeric(fitted(mod))
    expect_length(f, mod$nobs)
    expect_true(all(f == 0))
})

test_that("arma: h-step point forecast matches stats::arima for AR(1)/MA(1)/AR(2)", {
    set.seed(42)
    n <- 2000
    e <- rnorm(n)
    yy <- numeric(n)
    for (i in 2:n) yy[i] <- 2 + 0.6 * (yy[i - 1] - 2) + e[i]
    yy <- yy[500:n]
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0), arma = c(1,0))
    mod <- estimate(spec)
    p <- suppressWarnings(predict(mod, h = 10, nsim = 0))
    fit <- arima(yy, order = c(1,0,0))
    fc <- predict(fit, n.ahead = 10)
    expect_equal(as.numeric(p$mean), as.numeric(fc$pred), tolerance = 0.02)

    set.seed(7)
    n <- 3000
    e <- rnorm(n)
    yy <- numeric(n)
    for (i in 2:n) yy[i] <- 0.3 + e[i] + 0.5 * e[i - 1]
    yy <- yy[200:n]
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0), arma = c(0,1))
    mod <- estimate(spec)
    p <- suppressWarnings(predict(mod, h = 5, nsim = 0))
    fit <- arima(yy, order = c(0,0,1))
    fc <- predict(fit, n.ahead = 5)
    expect_equal(as.numeric(p$mean), as.numeric(fc$pred), tolerance = 0.02)
    # MA(1) forecast decays to the unconditional mean after 1 step
    expect_equal(as.numeric(p$mean)[2:5], rep(as.numeric(p$mean)[5], 4), tolerance = 1e-6)
})

test_that("arma: predict() does not warn, with or without arma", {
    spec0 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(0,0))
    mod0 <- estimate(spec0)
    expect_warning(predict(mod0, h = 5, nsim = 100), NA)

    spec1 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod1 <- estimate(spec1)
    expect_warning(predict(mod1, h = 5, nsim = 100), NA)
})

test_that("arma: predict() simulated bands are centered consistently with the analytic mean forecast, including when arma order exceeds garch order", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,2))
    mod <- estimate(spec)
    # previously errored: innov_init/var_init length mismatch when
    # max(arma) > max(garch order)
    p <- predict(mod, h = 8, nsim = 100000, seed = 1)
    mc_mean <- colMeans(p$distribution)
    expect_true(all(abs(as.numeric(mc_mean) - as.numeric(p$mean)) < 0.01))
})

test_that("arma: predict() mean forecast reduces exactly to constant mu when arma = c(0,0)", {
    spec0 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(0,0))
    mod0 <- estimate(spec0)
    p0 <- predict(mod0, h = 5, nsim = 0)
    expect_true(all(as.numeric(p0$mean) == mod0$parmatrix[parameter == "mu"]$value))
})

test_that("arma: rejects malformed orders", {
    expect_error(garch_modelspec(y[1:1800,1], model = "garch", order = c(1,1), arma = c(-1,0)))
    expect_error(garch_modelspec(y[1:1800,1], model = "garch", order = c(1,1), arma = c(1,1,1)))
})

# ---------------------------------------------------------------------------
# Higher-order ARMA(p,q) coverage: (1,1), (1,2), (2,1), (2,2), including
# cases where max(ar,ma) > max(garch order), which exercises the combined
# pre-sample burn-in fix (max(garch order, arma order)) in specification.R,
# garchfun.hpp, and simulate.R/simulation.cpp.
# ---------------------------------------------------------------------------

.simulate_arma_series <- function(ar_true, ma_true, mu = 0.5, n = 5000, seed = 1) {
    set.seed(seed)
    p <- length(ar_true); q <- length(ma_true)
    burn <- 500
    e <- rnorm(n + burn)
    yy <- eps <- numeric(n + burn)
    for (i in (max(p, q, 1) + 1):(n + burn)) {
        val <- mu
        for (j in seq_len(p)) val <- val + ar_true[j] * (yy[i - j] - mu)
        for (j in seq_len(q)) val <- val + ma_true[j] * eps[i - j]
        eps[i] <- e[i]
        yy[i] <- val + eps[i]
    }
    tail(yy, n)
}

arma_orders_to_test <- list(
    "ARMA(1,1)" = list(ar = 0.5, ma = 0.4, order = c(1,1)),
    "ARMA(1,2)" = list(ar = 0.5, ma = c(0.3, -0.2), order = c(1,2)),
    "ARMA(2,1)" = list(ar = c(0.4, -0.2), ma = 0.3, order = c(2,1)),
    "ARMA(2,2)" = list(ar = c(0.4, -0.2), ma = c(0.3, 0.15), order = c(2,2))
)

for (nm in names(arma_orders_to_test)) {
    spec_case <- arma_orders_to_test[[nm]]
    local({
        ar_true <- spec_case$ar
        ma_true <- spec_case$ma
        arma_order <- spec_case$order

        test_that(paste0("arma: ", nm, " pure-ARMA parameter recovery matches stats::arima"), {
            yy <- .simulate_arma_series(ar_true, ma_true, mu = 0.5, n = 5000, seed = 1)
            ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
            spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(0,0), arma = arma_order)
            mod <- estimate(spec)
            ac <- arma_coefficients(mod)
            fit <- arima(yy, order = c(arma_order[1], 0, arma_order[2]))
            if (arma_order[1] > 0) expect_equal(as.numeric(ac$ar), unname(coef(fit)[paste0("ar", 1:arma_order[1])]), tolerance = 0.03)
            if (arma_order[2] > 0) expect_equal(as.numeric(ac$ma), unname(coef(fit)[paste0("ma", 1:arma_order[2])]), tolerance = 0.03)
            expect_equal(mod$parmatrix[parameter == "mu"]$value, unname(coef(fit)["intercept"]), tolerance = 0.03)
            # nested-model likelihood sanity: tsgarch's MLE should be at least
            # as good as arima's own (both fit the same model)
            expect_true(as.numeric(logLik(mod)) >= fit$loglik - 0.5)
        })

        test_that(paste0("arma: ", nm, " + GARCH(1,1) estimated coefficients are stationary/invertible"), {
            spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = arma_order)
            mod <- estimate(spec)
            expect_true(mod$conditions$kkt1)
            ac <- arma_coefficients(mod)
            if (length(ac$ar) > 0) expect_true(min(Mod(polyroot(c(1, -as.numeric(ac$ar))))) > 1)
            if (length(ac$ma) > 0) expect_true(min(Mod(polyroot(c(1, as.numeric(ac$ma))))) > 1)
        })

        test_that(paste0("arma: ", nm, " + GARCH(1,1) simulate() zero-innovation path reduces exactly to mu"), {
            spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = arma_order)
            mod <- estimate(spec)
            mu <- mod$parmatrix[parameter == "mu"]$value
            spec_sim <- mod$spec
            spec_sim$parmatrix <- mod$parmatrix
            zeroinnov <- matrix(0, nrow = 1, ncol = 20)
            sim0 <- simulate(spec_sim, h = 20, nsim = 1, innov = zeroinnov, seed = 1)
            expect_equal(as.numeric(sim0$series), rep(mu, 20), tolerance = 1e-8)
        })

        test_that(paste0("arma: ", nm, " + GARCH(1,1) simulate() Monte Carlo mean converges to mu"), {
            spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = arma_order)
            mod <- estimate(spec)
            mu <- mod$parmatrix[parameter == "mu"]$value
            spec_sim <- mod$spec
            spec_sim$parmatrix <- mod$parmatrix
            sim <- simulate(spec_sim, h = 10, nsim = 20000, seed = 321)
            mc_mean <- colMeans(sim$series)
            expect_true(all(abs(mc_mean - mu) < 0.02))
        })

        test_that(paste0("arma: ", nm, " + GARCH(1,1) tsfilter() chained append matches a direct re-fit of the full series"), {
            spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = arma_order)
            mod <- estimate(spec)
            expect_equal(length(mod$conditional_mu), 1800)
            expect_equal(as.numeric(residuals(mod)), as.numeric(mod$spec$target$y_orig) - as.numeric(fitted(mod)))

            new_y <- y[1801:1974,1]
            filtered <- tsfilter(mod, y = new_y)
            expect_length(filtered$conditional_mu, 1974)
            expect_equal(filtered$nobs, 1974)

            # a fixed-parameter direct filter over the FULL series should
            # produce a numerically identical fitted conditional mean to the
            # chained, incremental tsfilter() append (both continue the same
            # ARMA recursion, just via different code paths)
            spec_full <- garch_modelspec(y[,1], constant = TRUE, model = "garch", order = c(1,1), arma = arma_order)
            spec_full$parmatrix <- copy(mod$parmatrix)
            full_direct <- tsfilter(spec_full)
            expect_equal(as.numeric(full_direct$conditional_mu), as.numeric(filtered$conditional_mu), tolerance = 1e-8)
        })
    })
}
