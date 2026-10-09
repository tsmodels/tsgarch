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
