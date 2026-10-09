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

test_that(".parametric_standardized_residual_draws returns a valid n x B matrix", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod <- estimate(spec)
    z <- tsgarch:::.parametric_standardized_residual_draws(mod, B = 15)
    expect_equal(dim(z), c(mod$nobs, 15))
    expect_false(any(is.na(z)))
    # different columns (different draws) should not be identical
    expect_false(isTRUE(all.equal(z[,1], z[,2])))
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

test_that("plot arma panel errors informatively for unsupported envelope values", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod <- estimate(spec)
    expect_error(plot(mod, type = "arma", which = 3, envelope = "bogus"), "should be one of")
})

test_that("plot arma panel errors informatively for invalid which values", {
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch", order = c(1,1), arma = c(1,1))
    mod <- estimate(spec)
    expect_error(plot(mod, type = "arma", which = 5), "which must be an integer vector")
    expect_error(plot(mod, type = "arma", which = 0), "which must be an integer vector")
})
