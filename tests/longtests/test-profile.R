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
