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

test_that("arma_near_cancellation detects close AR/MA roots", {
    roots <- list(ar = 0.5 + 0i, ma = 0.55 + 0i)
    out <- arma_near_cancellation(roots, tol = 0.1)
    expect_gt(nrow(out), 0)
    expect_lt(out$distance[1], 0.1)
})
