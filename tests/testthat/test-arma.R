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

test_that("arma: rejects malformed orders", {
    expect_error(garch_modelspec(y[1:1800,1], model = "garch", order = c(1,1), arma = c(-1,0)))
    expect_error(garch_modelspec(y[1:1800,1], model = "garch", order = c(1,1), arma = c(1,1,1)))
})

test_that("arma: ARMA(2,2) + GARCH(1,1) simulate() zero-innovation path reduces exactly to mu", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,2))
    mod <- estimate(spec)
    mu <- mod$parmatrix[parameter == "mu"]$value
    spec_sim <- mod$spec
    spec_sim$parmatrix <- mod$parmatrix
    zeroinnov <- matrix(0, nrow = 1, ncol = 20)
    sim0 <- simulate(spec_sim, h = 20, nsim = 1, innov = zeroinnov, seed = 1)
    expect_equal(as.numeric(sim0$series), rep(mu, 20), tolerance = 1e-8)
})

test_that("arma: ARMA(2,2) + GARCH(1,1) tsfilter() chained append matches a direct re-fit of the full series", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,2))
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
    spec_full <- garch_modelspec(y[,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,2))
    spec_full$parmatrix <- copy(mod$parmatrix)
    full_direct <- tsfilter(spec_full)
    expect_equal(as.numeric(full_direct$conditional_mu), as.numeric(filtered$conditional_mu), tolerance = 1e-8)
})

test_that("arma: summary()/print()/as_flextable() show the transformed ar/ma coefficients, not raw arpacf/mapacf", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,1))
    mod <- estimate(spec)
    s <- summary(mod)
    expect_true(all(c("ar1","ar2","ma1") %in% s$coefficients$term))
    expect_false(any(grepl("^arpacf|^mapacf", s$coefficients$term)))
    ac <- arma_coefficients(mod)
    expect_equal(s$coefficients[term == "ar1"]$Estimate, unname(ac$ar["ar1"]))
    expect_equal(s$coefficients[term == "ar2"]$Estimate, unname(ac$ar["ar2"]))
    expect_equal(s$coefficients[term == "ma1"]$Estimate, unname(ac$ma["ma1"]))
    expect_equal(s$symbol[s$coefficients$term == "ar1"], "\\phi_1")
    expect_equal(s$symbol[s$coefficients$term == "ma1"], "\\theta_1")
    # unaffected when arma = c(0,0)
    spec0 <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1))
    mod0 <- estimate(spec0)
    s0 <- summary(mod0)
    expect_false(any(grepl("^ar[0-9]|^ma[0-9]", s0$coefficients$term)))
})

test_that("arma: is allowed for all 8 native variants including igarch/ewma", {
    for (mdl in c("garch", "egarch", "gjrgarch", "aparch", "fgarch", "cgarch", "igarch", "ewma")) {
        expect_error(garch_modelspec(y[1:1800,1], constant = TRUE, model = mdl, order = c(1,1), arma = c(1,1)), NA)
    }
})

# ---------------------------------------------------------------------------
# Fixing an entire AR and/or MA polynomial at target coefficients (set
# value = target coefficients, estimate = 0, on every row of the group).
# ---------------------------------------------------------------------------

test_that("arma: fixing the entire ar polynomial (order 3) recovers the exact joint target despite lag coupling", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(3,1))
    target_ar <- c(0.3, -0.1, 0.05)
    spec$parmatrix[group == "arpacf", value := target_ar]
    spec$parmatrix[group == "arpacf", estimate := 0]
    mod <- estimate(spec)
    ac <- arma_coefficients(mod)
    expect_equal(unname(ac$ar), target_ar, tolerance = 1e-6)
    # ma (order 1, no coupling with ar) remains free and estimated
    expect_true("ma1" %in% summary(mod)$coefficients$term)
})

test_that("arma: partially fixing only some lags of a polynomial is rejected", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(3,1))
    spec$parmatrix[parameter == "arpacf1", estimate := 0]
    expect_error(estimate(spec), "partial fixing")
})

test_that("arma: fixing a non-stationary AR target or non-invertible MA target is rejected", {
    spec_ar <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    spec_ar$parmatrix[parameter == "arpacf1", value := 1.5]
    spec_ar$parmatrix[parameter == "arpacf1", estimate := 0]
    expect_error(estimate(spec_ar), "stationary")

    spec_ma <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    spec_ma$parmatrix[parameter == "mapacf1", value := -1.2]
    spec_ma$parmatrix[parameter == "mapacf1", estimate := 0]
    expect_error(estimate(spec_ma), "invertible")
})

test_that("arma: bootstrap predictive distribution is centered on the ARMA mean forecast, not mu", {
    # persistent AR(1) mean + garch(1,1) errors; the last observation is
    # placed far from mu so the 1-step mean forecast and mu are far apart
    # (otherwise the test cannot discriminate the two centerings)
    set.seed(99)
    n <- 1500
    mu_true <- 0.5
    phi <- 0.8
    sigma2 <- eps <- yy <- numeric(n)
    sigma2[1] <- 0.2
    for (i in 2:n) {
        sigma2[i] <- 0.05 + 0.1 * eps[i - 1]^2 + 0.85 * sigma2[i - 1]
        eps[i] <- rnorm(1) * sqrt(sigma2[i])
        yy[i] <- mu_true + phi * (yy[i - 1] - mu_true) + eps[i]
    }
    yy[n] <- mu_true + 4
    ys <- xts(yy, as.Date(seq_along(yy), origin = "1970-01-01"))
    spec <- garch_modelspec(ys, model = "garch", constant = TRUE, order = c(1,1), arma = c(1,0))
    mod <- estimate(spec)
    mu <- mod$parmatrix[parameter == "mu"]$value
    for (m in c("parametric", "bootstrap")) {
        p <- predict(mod, h = 5, nsim = 3000, sim_method = m, seed = 42)
        mc_mean1 <- mean(p$distribution[,1])
        # sanity: the 1-step mean forecast must be far from mu
        expect_true(abs(as.numeric(p$mean)[1] - mu) > 1)
        expect_true(abs(mc_mean1 - as.numeric(p$mean)[1]) < 0.1)
        expect_true(abs(mc_mean1 - as.numeric(p$mean)[1]) < abs(mc_mean1 - mu))
    }
})

test_that("arma: variance_targeting constant_variance uses ARMA innovations eps^2 (estimate and filter agree with TMB)", {
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1),
                            arma = c(1,1), variance_targeting = TRUE)
    mod <- estimate(spec)
    expect_equal(as.numeric(unconditional(mod)), mean(as.numeric(residuals(mod))^2), tolerance = 1e-6)
    expect_equal(as.numeric(mod$target_omega), as.numeric(unconditional(mod)) * (1 - as.numeric(persistence(mod))), tolerance = 1e-6)

    # same consistency after a tsfilter() round-trip on the spec (filter.R copy)
    spec_f <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1),
                              arma = c(1,1), variance_targeting = TRUE)
    spec_f$parmatrix <- copy(mod$parmatrix)
    f <- tsfilter(spec_f)
    expect_equal(f$constant_variance, mean((as.numeric(f$spec$target$y_orig) - as.numeric(f$conditional_mu))^2), tolerance = 1e-8)
    expect_equal(as.numeric(unconditional(f)), mean(as.numeric(residuals(f))^2), tolerance = 1e-6)
    expect_equal(as.numeric(f$target_omega), as.numeric(unconditional(f)) * (1 - as.numeric(persistence(f))), tolerance = 1e-6)
})
