test_that("simulate: a non-positive implied initial variance errors rather than returning zeros",{
    # ewma fixes omega at zero, so initv = omega/(1 - 0.999) = 0 and the
    # variance recursion would stay at zero for every step, silently returning
    # sigma == 0 and a series equal to mu. var_init must be supplied instead.
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "ewma")
    mod <- estimate(spec)
    spec$parmatrix <- copy(mod$parmatrix)
    expect_error(simulate(spec, h = 5, nsim = 2, seed = 1), "var_init")
    sim <- simulate(spec, h = 5, nsim = 2, var_init = tail(as.numeric(sigma(mod)), 1)^2, seed = 1)
    expect_true(all(as.matrix(sim$sigma) > 0))
    # every other flavour has a positive implied seed and is unaffected
    for (m in c("garch","igarch","egarch","aparch","fgarch","gjrgarch","cgarch")) {
        sp <- garch_modelspec(y[1:1800,1], constant = TRUE, model = m)
        mo <- suppressWarnings(estimate(sp))
        sp$parmatrix <- copy(mo$parmatrix)
        s <- suppressWarnings(simulate(sp, h = 5, nsim = 2, seed = 1))
        expect_true(all(as.matrix(s$sigma) > 0), info = m)
    }
})

test_that("gjrgarch(2,3) simulation: validate algoritm",{
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch",
                                  order = c(2,3), arma = c(0,0), vreg = y[1:1800,2],
                                  distribution = "norm")
    local_mod <- suppressWarnings(estimate(local_spec))
    spec <- copy(local_spec)
    spec$parmatrix <- copy(local_mod$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod, standardize = TRUE)), nrow = 1)
    maxpq <- max(local_spec$model$order, local_spec$model$arma)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod$var_initial,
                    innov = z, vreg = v, innov_init = rep(1, maxpq),
                    arch_initial = local_mod$arch_initial)
    expect_equal(sim$sigma[1,], local_mod$sigma)
})

test_that("aparch(1,1) arma(2,1) simulation: validate algoritm",{
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                  order = c(1,1), arma = c(2,1), vreg = y[1:1800,2],
                                  distribution = "norm")
    local_mod <- suppressWarnings(estimate(local_spec))
    spec <- copy(local_spec)
    spec$parmatrix <- copy(local_mod$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod, standardize = TRUE)), nrow = 1)
    maxpq <- max(local_spec$model$order, local_spec$model$arma)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod$var_initial,
                    innov = z, vreg = v, innov_init = rep(1, maxpq),
                    arch_initial = local_mod$arch_initial)
    expect_equal(sim$sigma[1,], local_mod$sigma)
})

test_that("aparch(2,3) simulation: validate algoritm",{
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                  order = c(2,3), arma = c(0,0), vreg = y[1:1800,2],
                                  distribution = "norm")
    local_mod <- suppressWarnings(estimate(local_spec))
    spec <- copy(local_spec)
    spec$parmatrix <- copy(local_mod$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod, standardize = TRUE)), nrow = 1)
    maxpq <- max(local_spec$model$order, local_spec$model$arma)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod$var_initial,
                    innov = z, vreg = v, innov_init = rep(1, maxpq),
                    arch_initial = local_mod$arch_initial)
    expect_equal(sim$sigma[1,], local_mod$sigma)
})

test_that("fgarch(1,1) arma(2,1) simulation: validate algoritm",{
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch",
                                  order = c(1,1), arma = c(2,1), vreg = y[1:1800,2],
                                  distribution = "norm")
    local_mod <- suppressWarnings(estimate(local_spec))
    spec <- copy(local_spec)
    spec$parmatrix <- copy(local_mod$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod, standardize = TRUE)), nrow = 1)
    maxpq <- max(local_spec$model$order, local_spec$model$arma)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod$var_initial,
                    innov = z, vreg = v, innov_init = rep(1, maxpq),
                    arch_initial = local_mod$arch_initial)
    expect_equal(sim$sigma[1,], local_mod$sigma)
})

test_that("fgarch(2,3) simulation: validate algoritm",{
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch",
                                  order = c(2,3), arma = c(0,0), vreg = y[1:1800,2],
                                  distribution = "norm")
    local_mod <- suppressWarnings(estimate(local_spec))
    spec <- copy(local_spec)
    spec$parmatrix <- copy(local_mod$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod, standardize = TRUE)), nrow = 1)
    maxpq <- max(local_spec$model$order, local_spec$model$arma)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod$var_initial,
                    innov = z, vreg = v, innov_init = rep(1, maxpq),
                    arch_initial = local_mod$arch_initial)
    expect_equal(sim$sigma[1,], local_mod$sigma)
})

test_that("predict: egarch simulation branch works when arma order exceeds garch order",{
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "egarch",
                            order = c(2,1), arma = c(3,0), distribution = "norm")
    mod <- suppressWarnings(estimate(spec))
    p <- predict(mod, h = 5, nsim = 100, seed = 1)
    expect_length(as.numeric(p$sigma), 5)
    expect_true(all(is.finite(as.numeric(p$sigma)) & as.numeric(p$sigma) > 0))
})
