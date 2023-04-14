
test_that("simulate norm: same seed same output",{
    spec <- global_spec_norm
    spec$parmatrix <- copy(global_mod_norm$parmatrix)
    maxpq <- max(spec$model$order)
    v_init <- as.numeric(tail(sigma(global_mod_norm)^2, maxpq))
    i_init <- as.numeric(tail(residuals(global_mod_norm),maxpq))
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
    spec <- global_spec_ghst
    spec$parmatrix <- copy(global_mod_ghst$parmatrix)
    maxpq <- max(spec$model$order)
    v_init <- as.numeric(tail(sigma(global_mod_ghst)^2, maxpq))
    i_init <- as.numeric(tail(residuals(global_mod_ghst),maxpq))
    simulate_spec1 <- simulate(spec, nsim = 100, seed = 101, h = 10, var_init = v_init,
                               innov_init = i_init, vreg = y[1801:1810,2])
    simulate_spec2 <- simulate(spec, nsim = 100, seed = 101, h = 10, var_init = v_init,
                               innov_init = i_init, vreg = y[1801:1810,2])
    expect_equal(simulate_spec1$series,simulate_spec2$series)
})
