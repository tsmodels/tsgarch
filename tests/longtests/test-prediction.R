test_that("prediction long-run: GARCH(1,1)",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                         order = c(1,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
})

test_that("prediction long-run: GARCH(2,1)",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                        order = c(2,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
})

test_that("prediction long-run: EGARCH(1,1)-NORM",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch",
                                        order = c(1,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 1000, seed = 100, burn = 1000)
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u, tolerance = 0.0001)
    expect_equal(mean(rowMeans(local_mod_sim$sigma^2), trim = 0.01), u, tolerance = 0.01)
})

test_that("prediction long-run: EGARCH(1,1)-GHST",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch",
                                        order = c(1,1), distribution = "ghst")
    local_mod_garch <- suppressWarnings(estimate(local_spec_garch))
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 1000, seed = 100, burn = 1000)
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u, tolerance = 0.0001)
    expect_equal(mean(rowMeans(local_mod_sim$sigma^2), trim = 0.01), u, tolerance = 0.01)
})

test_that("prediction long-run: EGARCH(1,1)-GH",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch",
                                        order = c(1,1), distribution = "gh")
    local_mod_garch <- estimate(local_spec_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 1000, seed = 100, burn = 1000)
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u, tolerance = 0.0001)
    expect_equal(mean(rowMeans(local_mod_sim$sigma^2), trim = 0.01), u, tolerance = 0.01)
})

test_that("prediction long-run: EGARCH(2,1)-NORM",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch",
                                        order = c(2,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 1000, seed = 100, burn = 1000)
    p <- predict(local_mod_garch, h = 4000, nsim = 0, seed = 500)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u, tolerance = 0.01)
    expect_equal(mean(rowMeans(local_mod_sim$sigma^2)), u, tolerance = 0.01)
})


test_that("prediction long-run: APARCH(1,1)-NORM",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                        order = c(1,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    delta <- coef(local_mod_garch)["delta"]
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 5000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^delta))^(2/delta))
    expect_equal(sim_u, u, tolerance = 0.001)
})


test_that("prediction long-run: APARCH(1,1)-GH",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                        order = c(1,1), distribution = "gh")
    local_mod_garch <- estimate(local_spec_garch)
    delta <- coef(local_mod_garch)["delta"]
    p <- predict(local_mod_garch, h = 1500)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 5000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^delta))^(2/delta))
    expect_equal(sim_u, u, tolerance = 0.01)
})

test_that("prediction long-run: APARCH(1,1)-JSU",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                        order = c(1,1), distribution = "jsu")
    local_mod_garch <- estimate(local_spec_garch)
    delta <- coef(local_mod_garch)["delta"]
    p <- predict(local_mod_garch, h = 1500)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 5000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^delta))^(2/delta))
    expect_equal(sim_u, u, tolerance = 0.01)
})


test_that("prediction long-run: APARCH(2,1)",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                        order = c(2,1), distribution = "norm")
    local_spec_garch$parmatrix[parameter == "alpha2", lower := -0.2]
    local_spec_garch$parmatrix[parameter == "gamma2", lower := -0.2]
    local_mod_garch <- estimate(local_spec_garch)
    delta <- coef(local_mod_garch)["delta"]
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 5000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u, tolerance = 0.001)
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^delta))^(2/delta))
    expect_equal(sim_u, u, tolerance = 0.01)
})

test_that("prediction long-run: CGARCH(1,1)",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch",
                                        order = c(1,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    p <- predict(local_mod_garch, h = 1500)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 3000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), unname(u), tolerance = 0.0001)
    expect_equal(mean(rowMeans(local_mod_sim$sigma^2)), unname(u), tolerance = 0.01)
})

test_that("prediction long-run: CGARCH(2,1)",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch",
                                        order = c(2,1))
    local_spec_garch$parmatrix[parameter == "alpha1", lower := -0.5]
    local_spec_garch$parmatrix[parameter == "alpha2", lower := -0.5]
    local_mod_garch <- estimate(local_spec_garch)
    p <- predict(local_mod_garch, h = 1500)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 3000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), unname(u), tolerance = 0.0001)
    expect_equal(mean(rowMeans(local_mod_sim$sigma^2)), unname(u), tolerance = 0.01)
})

test_that("prediction long-run: FGARCH(1,1)-NORM",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch",
                                        order = c(1,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    delta <- coef(local_mod_garch)["delta"]
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 5000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
    delta <- coef(local_mod_garch)["delta"]
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^delta))^(2/delta))
    expect_equal(sim_u, u, tolerance = 0.01)
})


test_that("prediction long-run: FGARCH(1,1)-SSTD",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch",
                                        order = c(1,1), distribution = "sstd")
    local_mod_garch <- estimate(local_spec_garch)
    p <- predict(local_mod_garch, h = 1500)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 8000, seed = 800, burn = 2000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
    delta <- coef(local_mod_garch)["delta"]
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^delta))^(2/delta))
    expect_equal(sim_u, u, tolerance = 0.03)
})


test_that("prediction long-run: GJR(1,1)-NORM",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch",
                                        order = c(1,1), distribution = "norm")
    local_mod_garch <- estimate(local_spec_garch)
    p <- predict(local_mod_garch, h = 1000)
    u <- unconditional(local_mod_garch)
    spec_copy <- local_spec_garch
    spec_copy$parmatrix <- copy(local_mod_garch$parmatrix)
    local_mod_sim <- simulate(spec_copy, nsim = 1000, h = 5000, seed = 100, burn = 1000)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u)
    sim_u <- unname(mean(rowMeans(local_mod_sim$sigma^2)))
    expect_equal(sim_u, u, tolerance = 0.001)
})

test_that("prediction long-run: GJR(1,1)-SSTD",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch",
                                        order = c(1,1), distribution = "sstd")
    local_mod_garch <- estimate(local_spec_garch, stationarity_constraint = 0.97)
    p <- predict(local_mod_garch, h = 4500)
    u <- unconditional(local_mod_garch)
    expect_equal(as.numeric(tail(p$sigma^2,1)), u, tolerance = 0.00001)
})
