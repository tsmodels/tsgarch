test_that("garch(1,1) profile",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                         order = c(1,1), distribution = "jsu")
    mod <- estimate(spec)
    new_spec <- copy(spec)
    new_spec$parmatrix <- copy(mod$parmatrix)
    p <- suppressWarnings(tsprofile(new_spec, nsim = 10, sizes = c(500, 1000), var_init = mod$var_initial, seed = 100))
    rmse_500 <- mean(p$summary[size == 500]$RMSE)
    rmse_1000 <- mean(p$summary[size == 1000]$RMSE)
    expect_equal(rmse_500, 0.2766035, tolerance =  0.0001)
    expect_equal(rmse_1000, 0.1790083, tolerance = 0.0001)
    expect_length(colnames(p$summary), 9)
})

test_that("profile propagates the arma order and refuses mean regressors",{
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch",
                            order = c(1,1), arma = c(1,0))
    mod <- estimate(spec)
    new_spec <- copy(spec)
    new_spec$parmatrix <- copy(mod$parmatrix)
    p <- suppressWarnings(tsprofile(new_spec, nsim = 2, sizes = 400, seed = 100))
    # before the fix the re-fitted specs were constant-mean, so the ARMA
    # parameter was never estimated and never appeared in the profile
    expect_true("arpacf1" %in% as.character(p$summary$parameter))
    expect_true(is.finite(p$summary[parameter == "arpacf1"]$MEAN[1]))
    # mean equation regressors are refused rather than silently profiled with
    # a zero regressor contribution in the simulated data
    xspec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "garch",
                             order = c(1,1), xreg = y[1:800,2])
    xmod <- estimate(xspec)
    new_xspec <- copy(xspec)
    new_xspec$parmatrix <- copy(xmod$parmatrix)
    expect_error(tsprofile(new_xspec, nsim = 2, sizes = 400), "mean equation")
})
