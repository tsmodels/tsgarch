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
                    var_init = rep(global_mod_cgarch$var_initial,2),
                    innov = z, vreg = v)
    expect_equal(sim$sigma[1,], local_mod_cgarch$sigma)

})

test_that("egarch(1,1) simulation: validate algoritm",{
    spec <- copy(global_spec_garch)
    spec$parmatrix <- copy(global_mod_garch$parmatrix)
    v <- c(as.numeric(y[1:1800,2]) * coef(global_mod_garch)["xi1"])
    z <- matrix(as.numeric(residuals(global_mod_garch, standardize = TRUE)), nrow = 1)
    # use fixed innovation and replicate the initial conditions to guarantee a deterministic
    # simulation which serves to validate the algorithm for correctness and reproducability
    sim <- simulate(spec, nsim = 1, h = length(spec$target$y_orig),
                    var_init = global_mod_garch$var_initial,
                    innov = z, vreg = v,
                    arch_initial = global_mod_garch$arch_initial)
    expect_equal(sim$sigma[1,], global_mod_garch$sigma)
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
    expect_equal(simulate_spec1$series,simulate_spec2$series)
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
    expect_equal(simulate_spec1$series,simulate_spec2$series)
})

test_that("simulation: long run variance check",{
    spec_garch <- garch_modelspec(y = y[1:1800,1], constant = TRUE, model = "garch")
    sim <- simulate(spec_garch, nsim = 1, h = 25000, seed = 77, burn = 100)
    expect_equal(mean(sim$sigma[1,]^2), unconditional(spec_garch), tolerance = 0.01)
    spec_egarch <- garch_modelspec(y = y[1:1800,1], constant = TRUE, model = "egarch")
    sim <- simulate(spec_egarch, nsim = 1, h = 25000, seed = 727, burn = 100)
    expect_equal(mean(sim$sigma[1,]^2), unconditional(spec_egarch), tolerance = 0.1)

})
