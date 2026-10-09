test_that("garch(1,1) backtest:rolling",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), vreg = y[1:1800,2], distribution = "norm")
    b <- tsbacktest(spec, start = 1000, end = 1200, h = 1, estimate_every = 100, rolling = TRUE, trace = FALSE)
    f_dates_diff <- as.integer(unique(diff(b$table$forecast_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_sum <- as.integer(sum(diff(b$table$estimation_date)))
    expect_equal(f_dates_diff, 1)
    expect_equal(e_dates_diff, 100)
    expect_equal(e_dates_sum, 200 - 100)
    expect_equal(NROW(b$table), 200)
})

test_that("garch(1,1) backtest:rolling multihorizon",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), vreg = y[1:1800,2], distribution = "norm")
    b <- tsbacktest(spec, start = 1000, end = 1200, h = 5, estimate_every = 100, rolling = TRUE, trace = FALSE)
    f_dates_diff <- as.integer(unique(b$table$forecast_date - b$table$filter_date))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_sum <- as.integer(sum(diff(b$table$estimation_date)))
    expect_equal(max(f_dates_diff), 5)
    expect_equal(min(f_dates_diff), 1)
    expect_equal(e_dates_diff, 100)
    expect_equal(e_dates_sum, 200 - 100)
    expect_equal(NROW(b$table), 5 * 200 - 5*2)
})

test_that("garch(1,1) backtest: non rolling",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), vreg = y[1:1800,2], distribution = "norm")
    b <- tsbacktest(spec, start = 1000, end = 1200, h = 100, estimate_every = 100, rolling = FALSE, trace = FALSE)
    f_dates_diff <- as.integer(unique(diff(b$table$forecast_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_sum <- as.integer(sum(diff(b$table$estimation_date)))
    expect_equal(max(f_dates_diff), 1)
    expect_equal(min(f_dates_diff), 1)
    expect_equal(e_dates_diff, 100)
    expect_equal(e_dates_sum, 200 - 100)
    expect_equal(NROW(b$table), 200)
})


test_that("garch(1,1) backtest: non rolling non overlapping",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1), vreg = y[1:1800,2], distribution = "norm")
    b <- tsbacktest(spec, start = 1000, end = 1200, h = 80, estimate_every = 100, rolling = FALSE, trace = FALSE)
    f_dates_diff <- as.integer(unique(diff(b$table$forecast_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_diff <- as.integer(max(diff(b$table$estimation_date)))
    e_dates_sum <- as.integer(sum(diff(b$table$estimation_date)))
    expect_equal(max(f_dates_diff), 21)
    expect_equal(min(f_dates_diff), 1)
    expect_equal(e_dates_diff, 100)
    expect_equal(e_dates_sum, 200 - 100)
    expect_equal(NROW(b$table), 200 - 40)
})

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
