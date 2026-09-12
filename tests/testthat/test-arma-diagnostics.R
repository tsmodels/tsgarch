test_that("arma_inverse_roots and arma_irf on an estimated ARMA(2,1)-GARCH model", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,1))
    mod <- estimate(spec)
    ac <- arma_coefficients(mod)
    roots <- arma_inverse_roots(mod)

    expect_length(roots$ar, 2)
    expect_length(roots$ma, 1)
    # AR polynomial: 1 - ar1 z - ar2 z^2
    ar_roots_expected <- 1 / polyroot(c(1, -as.numeric(ac$ar)))
    expect_equal(sort(Mod(roots$ar)), sort(Mod(ar_roots_expected)), tolerance = 1e-10)
    # MA polynomial: 1 + ma1 z -> inverse root is -ma1
    expect_equal(as.numeric(roots$ma), -as.numeric(ac$ma), tolerance = 1e-10)

    irf <- arma_irf(mod)
    # psi_0 = 1, psi decays, cumulative exists
    expect_equal(as.numeric(irf$psi)[1], 1)
    expect_length(irf$cumulative, length(irf$psi))
    expect_true(all(abs(irf$psi[-1]) < 1))
})

test_that("arma_irf for AR(1) follows phi^k", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,0))
    mod <- estimate(spec)
    ac <- arma_coefficients(mod)
    irf <- arma_irf(mod)
    expected <- (as.numeric(ac$ar) ^ (0:9))
    expect_equal(as.numeric(irf$psi)[1:10], expected, tolerance = 1e-3)
})

test_that("arma_irf for MA(1) is zero after lag 1", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(0,1))
    mod <- estimate(spec)
    ac <- arma_coefficients(mod)
    irf <- arma_irf(mod)
    expect_equal(as.numeric(irf$psi)[1:2], c(1, as.numeric(ac$ma)), tolerance = 1e-3)
    expect_true(all(abs(as.numeric(irf$psi)[-c(1,2)]) < 1e-3))
})

test_that("arma_inverse_roots returns zero-length vectors for a constant-mean model", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1))
    mod <- estimate(spec)
    roots <- arma_inverse_roots(mod)
    expect_length(roots$ar, 0)
    expect_length(roots$ma, 0)
})

test_that("arma_near_cancellation detects close AR/MA roots", {
    roots <- list(ar = 0.5 + 0i, ma = 0.55 + 0i)
    out <- arma_near_cancellation(roots, tol = 0.1)
    expect_gt(nrow(out), 0)
    expect_lt(out$distance[1], 0.1)
})

test_that("plot method dispatches to garch panel and preserves par state", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod <- estimate(spec)
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    tmp <- tempfile(fileext = ".png")
    png(tmp)
    out <- plot(mod)
    dev.off()
    expect_identical(out, mod)
    # verify par is restored to original state (with tolerance for unimportant fields)
    expect_equal(par(no.readonly = TRUE)$mfrow, old$mfrow)
    expect_equal(par(no.readonly = TRUE)$mar, old$mar)
})

test_that("plot method dispatches to arma panel and returns object", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(2,1))
    mod <- estimate(spec)
    tmp <- tempfile(fileext = ".png")
    png(tmp)
    out <- plot(mod, type = "arma")
    dev.off()
    expect_identical(out, mod)

    tmp2 <- tempfile(fileext = ".png")
    png(tmp2)
    out2 <- plot(mod, type = "arma", which = 1)
    dev.off()
    expect_identical(out2, mod)
})

test_that("plot type validation errors on invalid type", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod <- estimate(spec)
    expect_error(plot(mod, type = "bogus"), "should be one of")
})

test_that("plot type = 'arma' errors for models without ARMA", {
    spec0 <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1))
    mod0 <- estimate(spec0)
    expect_error(plot(mod0, type = "arma"), "ARMA mean equation")
})
