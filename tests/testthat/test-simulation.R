test_that("garch(1,1) simulation: validate algoritm",{
    spec <- copy(global_spec_garch)
    spec$parmatrix <- copy(global_mod_garch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(global_mod_garch)["xi1"])
    z <- matrix(as.numeric(residuals(global_mod_garch, standardize = TRUE)), nrow = 1)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim1 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = global_mod_garch$var_initial,
                    innov = z, vreg = v,
                    arch_initial = global_mod_garch$arch_initial)

    sim2 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = global_mod_garch$var_initial,
                    innov = z, vreg = v)

    expect_equal(sim1$sigma[1,], global_mod_garch$sigma)
    expect_equal(sim2$sigma[1,], global_mod_garch$sigma)
})


test_that("gjrgarch(1,1) simulation: validate algoritm",{
    spec <- copy(global_spec_gjrgarch)
    spec$parmatrix <- copy(global_mod_gjrgarch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(global_mod_gjrgarch)["xi1"])
    z <- matrix(as.numeric(residuals(global_mod_gjrgarch, standardize = TRUE)), nrow = 1)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim1 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_gjrgarch$var_initial,
                     innov = z, vreg = v, innov_init = 1,
                     arch_initial = global_mod_gjrgarch$arch_initial)

    sim2 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_gjrgarch$var_initial,
                     innov = z, vreg = v)

    expect_equal(sim1$sigma[1,], global_mod_gjrgarch$sigma)
    # after initial conditions die out, we converge
    # reason for not having the exact same outomce for sim2 is because
    # arch_initial is not passed (this is the mean of the squared residuals * sample_kappa)
    expect_equal(tail(sim2$sigma[1,],10), tail(global_mod_gjrgarch$sigma,10))
})

test_that("aparch(1,1) simulation: validate algoritm",{
    spec <- copy(global_spec_aparch)
    spec$parmatrix <- copy(global_mod_aparch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(global_mod_aparch)["xi1"])
    z <- matrix(as.numeric(residuals(global_mod_aparch, standardize = TRUE)), nrow = 1)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim1 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_aparch$var_initial,
                     innov = z, vreg = v, innov_init = 1,
                     arch_initial = global_mod_aparch$arch_initial)

    sim2 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_aparch$var_initial,
                     innov = z, vreg = v)

    sim3 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     innov = z, vreg = v)

    expect_equal(sim1$sigma[1,], global_mod_aparch$sigma)
    expect_equal(tail(sim2$sigma[1,],10), tail(global_mod_aparch$sigma,10))
    expect_equal(tail(sim3$sigma[1,],10), tail(global_mod_aparch$sigma,10))
})



test_that("fgarch(1,1) simulation: validate algoritm",{
    spec <- copy(global_spec_fgarch)
    spec$parmatrix <- copy(global_mod_fgarch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(global_mod_fgarch)["xi1"])
    z <- matrix(as.numeric(residuals(global_mod_fgarch, standardize = TRUE)), nrow = 1)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim1 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_fgarch$var_initial,
                     innov = z, vreg = v, innov_init = 1,
                     arch_initial = global_mod_fgarch$arch_initial)

    sim2 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_fgarch$var_initial,
                     innov = z, vreg = v)

    sim3 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     innov = z, vreg = v)

    expect_equal(sim1$sigma[1,], global_mod_fgarch$sigma)
    expect_equal(tail(sim2$sigma[1,],10), tail(global_mod_fgarch$sigma,10))
    expect_equal(tail(sim3$sigma[1,],10), tail(global_mod_fgarch$sigma,10))
})

test_that("garch(2,1) simulation: validate algoritm",{
    local_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                         order = c(2,1), vreg = y[1:1800,2],
                                         distribution = "norm")
    local_mod_garch <- suppressWarnings(estimate(local_spec_garch))
    spec <- copy(local_spec_garch)
    spec$parmatrix <- copy(local_mod_garch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod_garch)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod_garch, standardize = TRUE)), nrow = 1)
    sim1 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod_garch$var_initial,
                    innov = z, vreg = v,
                    arch_initial = local_mod_garch$arch_initial)
    sim2 <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = local_mod_garch$var_initial,
                     innov = z, vreg = v)

    expect_equal(sim1$sigma[1,], local_mod_garch$sigma)
    expect_equal(sim2$sigma[1,], local_mod_garch$sigma)

})


test_that("cgarch(1,1) simulation: validate algoritm",{
    spec <- copy(global_spec_cgarch)
    spec$parmatrix <- copy(global_mod_cgarch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(global_mod_cgarch)["xi1"])
    z <- matrix(as.numeric(residuals(global_mod_cgarch, standardize = TRUE)), nrow = 1)
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                     var_init = global_mod_cgarch$var_initial,
                     innov = z, vreg = v)
    expect_equal(sim$sigma[1,], global_mod_cgarch$sigma)
})

test_that("cgarch(1,1) simulation: validate algoritm",{
    local_spec_cgarch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch",
                                         order = c(1,1), vreg = y[1:1800,2], multiplicative = TRUE,
                                         distribution = "norm")
    local_mod_cgarch <- suppressWarnings(estimate(local_spec_cgarch))
    spec <- copy(local_spec_cgarch)
    spec$parmatrix <- copy(local_mod_cgarch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod_cgarch)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod_cgarch, standardize = TRUE)), nrow = 1)
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = local_mod_cgarch$var_initial,
                    innov = z, vreg = v)
    expect_equal(sim$sigma[1,], local_mod_cgarch$sigma)
})

test_that("cgarch(2,1) simulation: validate algoritm",{
    local_spec_cgarch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch",
                                         order = c(2,1), vreg = y[1:1800,2],
                                         distribution = "norm")
    local_mod_cgarch <- suppressWarnings(estimate(local_spec_cgarch))
    spec <- copy(local_spec_cgarch)
    spec$parmatrix <- copy(local_mod_cgarch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(local_mod_cgarch)["xi1"])
    z <- matrix(as.numeric(residuals(local_mod_cgarch, standardize = TRUE)), nrow = 1)
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = rep(local_mod_cgarch$var_initial,2),
                    innov = z, vreg = v)
    expect_equal(sim$sigma[1,], local_mod_cgarch$sigma, tolerance = 0.01)

})

test_that("egarch(1,1) simulation: validate algoritm",{
    # this copied the garch globals throughout, so egarch was never actually
    # validated here; it now builds and fits an egarch spec
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch",
                                  order = c(1,1), arma = c(0,0), vreg = y[1:1800,2],
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

test_that("simulate norm: same seed same output",{
    spec <- copy(global_spec_garch)
    spec$parmatrix <- copy(global_mod_garch$parmatrix)
    maxpq <- max(spec$model$order)
    v_init <- as.numeric(tail(sigma(global_mod_garch)^2, maxpq))
    i_init <- as.numeric(tail(residuals(global_mod_garch),maxpq))
    simulate_spec1 <- simulate(spec, nsim = 100, seed = 101, h = 10, var_init = v_init,
                               innov_init = i_init, vreg = y[1801:1810,2])
    simulate_spec2 <- simulate(spec, nsim = 100, seed = 101, h = 10, var_init = v_init,
                               innov_init = i_init, vreg = y[1801:1810,2])
    expect_equal(simulate_spec1$series,simulate_spec2$series, tolerance = 0.001)
    expect_equal(NROW(simulate_spec1$series),100)
    expect_equal(NCOL(simulate_spec1$series),10)
    expect_s3_class(simulate_spec1, class = "tsgarch.simulate")
    expect_s3_class(simulate_spec1$sigma, class = "tsmodel.distribution")
})

test_that("simulate ghst: same seed same output",{
    spec <- global_spec_garch_jsu
    spec$parmatrix <- copy(global_mod_garch_jsu$parmatrix)
    maxpq <- max(spec$model$order)
    v_init <- as.numeric(tail(sigma(global_mod_garch_jsu)^2, maxpq))
    i_init <- as.numeric(tail(residuals(global_mod_garch_jsu),maxpq))
    simulate_spec1 <- simulate(spec, nsim = 100, seed = 101, h = 10, var_init = v_init,
                               innov_init = i_init, vreg = y[1801:1810,2])
    simulate_spec2 <- simulate(spec, nsim = 100, seed = 101, h = 10, var_init = v_init,
                               innov_init = i_init, vreg = y[1801:1810,2])
    expect_equal(simulate_spec1$series,simulate_spec2$series, tolerance = 0.001)
})

test_that("simulation: long run variance check",{
    spec_garch <- garch_modelspec(y = y[1:1800,1], constant = TRUE, model = "garch")
    sim <- simulate(spec_garch, nsim = 1, h = 25000, seed = 77, burn = 100)
    expect_equal(mean(sim$sigma[1,]^2), unconditional(spec_garch), tolerance = 0.01)
    spec_egarch <- garch_modelspec(y = y[1:1800,1], constant = TRUE, model = "egarch")
    sim <- simulate(spec_egarch, nsim = 1, h = 25000, seed = 727, burn = 100)
    expect_equal(mean(sim$sigma[1,]^2), unconditional(spec_egarch), tolerance = 0.1)

})

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

test_that("gjrgarch(1,1) arma(2,1) simulation: validate algoritm",{
    local_spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch",
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

test_that("simulate: identical innovation rows give identical sigma rows across flavours",{
    set.seed(1); one <- rnorm(60)
    Z <- matrix(rep(one, 3), nrow = 3, byrow = TRUE)
    # asymmetric innov_init: the squaring/abs in several init formulas would
    # mask a column-major scatter of a symmetric vector
    flavours <- c("garch","egarch","gjrgarch","aparch","fgarch","cgarch","igarch")
    for (m in flavours) {
        spec <- suppressWarnings(garch_modelspec(y[1:500,1], constant = TRUE, model = m,
                                  order = c(1,1), arma = c(2,1), distribution = "norm"))
        maxpq <- max(spec$model$order, spec$model$arma)
        ii <- seq(0.6, by = -1.1, length.out = maxpq)
        # var_init is subject to the same broadcast contract as innov_init, and
        # it has to vary over the pre-sample for a column-major fill to show up.
        # cgarch takes a maxpq by 2 matrix because it carries two variances: the
        # permanent (long run) component in column 1, and the total conditional
        # variance in column 2
        vi <- seq(0.9, by = 0.7, length.out = maxpq)
        if (m == "cgarch") vi <- cbind(vi, vi + 0.3)
        s1 <- suppressWarnings(simulate(spec, nsim = 1, h = 60, innov = matrix(one, nrow = 1), innov_init = ii, var_init = vi)$sigma)
        s3 <- suppressWarnings(simulate(spec, nsim = 3, h = 60, innov = Z, innov_init = ii, var_init = vi)$sigma)
        expect_equal(s1[1,], s3[1,], info = m)
        expect_equal(s3[1,], s3[2,], info = m)
        expect_equal(s3[1,], s3[3,], info = m)
    }
    # the same invariants with unequal per-lag parameters expose per-lag
    # recycling of gamma/eta across the pre-sample columns
    for (m in flavours) {
        spec <- suppressWarnings(garch_modelspec(y[1:500,1], constant = TRUE, model = m,
                                  order = c(2,1), arma = c(0,0), distribution = "norm"))
        for (g in c("alpha","gamma","eta","beta")) {
            idx <- which(spec$parmatrix$group == g)
            if (length(idx) > 1) spec$parmatrix[idx, value := value * seq(0.6, 1.4, length.out = length(idx))]
        }
        maxpq <- max(spec$model$order, spec$model$arma)
        ii <- seq(0.6, by = -1.1, length.out = maxpq)
        vi <- seq(0.9, by = 0.7, length.out = maxpq)
        if (m == "cgarch") vi <- cbind(vi, vi + 0.3)
        s1 <- suppressWarnings(simulate(spec, nsim = 1, h = 60, innov = matrix(one, nrow = 1), innov_init = ii, var_init = vi)$sigma)
        s3 <- suppressWarnings(simulate(spec, nsim = 3, h = 60, innov = Z, innov_init = ii, var_init = vi)$sigma)
        expect_equal(s1[1,], s3[1,], info = m)
        expect_equal(s3[1,], s3[2,], info = m)
        expect_equal(s3[1,], s3[3,], info = m)
    }
})

test_that("simulate: egarch default arch initialization matches the likelihood",{
    # egarchfun.hpp zeroes initial_arch, so the pre-sample arch term contributes
    # exactly nothing; with no innov_init to evaluate the equation at, the
    # simulation default must agree rather than use |0| - kappa = -kappa
    set.seed(1); one <- rnorm(40)
    spec <- suppressWarnings(garch_modelspec(y[1:500,1], constant = TRUE, model = "egarch",
                              order = c(2,1), arma = c(0,0), distribution = "norm"))
    idx <- which(spec$parmatrix$group == "gamma")
    spec$parmatrix[idx, value := c(0.12, 0.2)]
    maxpq <- max(spec$model$order, spec$model$arma)
    sdef <- suppressWarnings(simulate(spec, nsim = 1, h = 40, innov = matrix(one, nrow = 1))$sigma)
    szero <- suppressWarnings(simulate(spec, nsim = 1, h = 40, innov = matrix(one, nrow = 1),
                                       arch_initial = rep(0, maxpq))$sigma)
    expect_equal(as.numeric(sdef), as.numeric(szero))
    # the previous default, to show the check above is not vacuous
    skappa <- suppressWarnings(simulate(spec, nsim = 1, h = 40, innov = matrix(one, nrow = 1),
                                        arch_initial = rep(-sqrt(2/pi), maxpq))$sigma)
    expect_false(isTRUE(all.equal(as.numeric(sdef), as.numeric(skappa))))
})

test_that("predict: egarch simulation branch works when arma order exceeds garch order",{
    spec <- garch_modelspec(y[1:800,1], constant = TRUE, model = "egarch",
                            order = c(2,1), arma = c(3,0), distribution = "norm")
    mod <- suppressWarnings(estimate(spec))
    p <- predict(mod, h = 5, nsim = 100, seed = 1)
    expect_length(as.numeric(p$sigma), 5)
    expect_true(all(is.finite(as.numeric(p$sigma)) & as.numeric(p$sigma) > 0))
})
