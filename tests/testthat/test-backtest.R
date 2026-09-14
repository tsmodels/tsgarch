test_that("backtest: propagates the arma mean equation to rolling refits", {
    yb <- y[1:350,1]
    spec <- garch_modelspec(yb, model = "garch", constant = TRUE, order = c(1,1),
                            arma = c(1,0), distribution = "norm")
    b <- tsbacktest(spec, start = 347, end = 350, h = 1, rolling = FALSE, trace = FALSE)
    # the first estimation window trains on yb[1:347]
    refit <- estimate(garch_modelspec(yb[1:347], model = "garch", constant = TRUE,
                                      order = c(1,1), arma = c(1,0), distribution = "norm"))
    refit0 <- estimate(garch_modelspec(yb[1:347], model = "garch", constant = TRUE,
                                       order = c(1,1), arma = c(0,0), distribution = "norm"))
    p <- predict(refit, h = 1, nsim = 0)
    p0 <- predict(refit0, h = 1, nsim = 0)
    first <- b$table[1]
    expect_equal(first$mu, as.numeric(p$mean), tolerance = 1e-6)
    expect_false(isTRUE(all.equal(first$mu, as.numeric(p0$mean), tolerance = 1e-6)))
})

test_that("backtest: propagates the garch model flavor to rolling refits", {
    yb <- y[1:350,1]
    spec <- garch_modelspec(yb, model = "egarch", constant = TRUE, order = c(1,1),
                            distribution = "norm")
    b <- tsbacktest(spec, start = 347, end = 350, h = 1, rolling = FALSE, trace = FALSE)
    refit <- estimate(garch_modelspec(yb[1:347], model = "egarch", constant = TRUE,
                                      order = c(1,1), distribution = "norm"))
    refit0 <- estimate(garch_modelspec(yb[1:347], model = "garch", constant = TRUE,
                                       order = c(1,1), distribution = "norm"))
    p <- predict(refit, h = 1, nsim = 0)
    p0 <- predict(refit0, h = 1, nsim = 0)
    first <- b$table[1]
    expect_equal(first$sigma, as.numeric(p$sigma), tolerance = 1e-6)
    expect_false(isTRUE(all.equal(first$sigma, as.numeric(p0$sigma), tolerance = 1e-6)))
})

test_that("backtest: ewma refits stay ewma (not igarch)", {
    yb <- y[1:350,1]
    spec <- garch_modelspec(yb, model = "ewma", constant = TRUE, order = c(1,1),
                            distribution = "norm")
    b <- tsbacktest(spec, start = 347, end = 350, h = 1, rolling = FALSE, trace = FALSE)
    refit <- estimate(garch_modelspec(yb[1:347], model = "ewma", constant = TRUE,
                                      order = c(1,1), distribution = "norm"))
    refit0 <- estimate(garch_modelspec(yb[1:347], model = "igarch", constant = TRUE,
                                       order = c(1,1), distribution = "norm"))
    p <- predict(refit, h = 1, nsim = 0)
    p0 <- predict(refit0, h = 1, nsim = 0)
    first <- b$table[1]
    expect_equal(first$sigma, as.numeric(p$sigma), tolerance = 1e-6)
    expect_false(isTRUE(all.equal(first$sigma, as.numeric(p0$sigma), tolerance = 1e-6)))
})
