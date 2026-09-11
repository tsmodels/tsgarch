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

test_that("arma: rejects malformed orders", {
    expect_error(garch_modelspec(y[1:1800,1], model = "garch", order = c(1,1), arma = c(-1,0)))
    expect_error(garch_modelspec(y[1:1800,1], model = "garch", order = c(1,1), arma = c(1,1,1)))
})
