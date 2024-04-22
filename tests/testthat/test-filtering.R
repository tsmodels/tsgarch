
test_that("filter from estimate (garch) append: check results",{
    filt <- tsfilter(global_mod_garch, y = y[1801,1], newvreg = y[1801,2])
    expect_equal(filt$sigma[1:1800], global_mod_garch$sigma)
    expect_equal(residuals(global_mod_garch, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})

test_that("filter from spec (garch) no append: check results",{
    spec <- copy(global_spec_garch)
    spec$parmatrix <- copy(global_mod_garch$parmatrix)
    filt <- tsfilter(spec)
    expect_equal(filt$sigma[1:1800], global_mod_garch$sigma)
    expect_equal(residuals(global_mod_garch, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1800)
    expect_length(filt$spec$target$y_orig, 1800)
    expect_length(as.numeric(filt$spec$target$y), 1800)
    expect_equal(filt$nobs, 1800)
    expect_equal(filt$loglik, global_mod_garch$loglik)
})

test_that("filter from estimate (cgarch) append: check results",{
    filt <- tsfilter(global_mod_cgarch, y = y[1801,1], newvreg = y[1801,2])
    expect_equal(filt$sigma[1:1800], global_mod_cgarch$sigma)
    expect_equal(residuals(global_mod_cgarch, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})

test_that("filter from spec (cgarch) no append: check results",{
    spec <- copy(global_spec_cgarch)
    spec$parmatrix <- copy(global_mod_cgarch$parmatrix)
    filt <- tsfilter(spec)
    expect_equal(filt$sigma[1:1800], global_mod_cgarch$sigma)
    expect_equal(filt$permanent_component[1:1800], global_mod_cgarch$permanent_component)
    expect_equal(filt$transitory_component[1:1800], global_mod_cgarch$transitory_component)
    expect_equal(residuals(global_mod_cgarch, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1800)
    expect_length(filt$spec$target$y_orig, 1800)
    expect_length(as.numeric(filt$spec$target$y), 1800)
    expect_equal(filt$nobs, 1800)
    expect_equal(filt$loglik, global_mod_cgarch$loglik)
})


test_that("filter from estimate (gjrgarch) append: check results",{
    filt <- tsfilter(global_mod_gjrgarch, y = y[1801,1], newvreg = y[1801,2])
    expect_equal(filt$sigma[1:1800], global_mod_gjrgarch$sigma)
    expect_equal(residuals(global_mod_gjrgarch, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})

test_that("filter from spec (gjrgarch) no append: check results",{
    spec <- copy(global_spec_cgarch)
    spec$parmatrix <- copy(global_mod_cgarch$parmatrix)
    filt <- tsfilter(spec)
    expect_equal(filt$sigma[1:1800], global_mod_cgarch$sigma)
    expect_equal(filt$permanent_component[1:1800], global_mod_cgarch$permanent_component)
    expect_equal(filt$transitory_component[1:1800], global_mod_cgarch$transitory_component)
    expect_equal(residuals(global_mod_cgarch, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1800)
    expect_length(filt$spec$target$y_orig, 1800)
    expect_length(as.numeric(filt$spec$target$y), 1800)
    expect_equal(filt$nobs, 1800)
    expect_equal(filt$loglik, global_mod_cgarch$loglik)
})


test_that("filter from spec (cgarch) append with initial conditions: check results",{
    # see whether using same sample initialization with different data sizes
    # provides the same results
    local_spec_cgarch <- garch_modelspec(y[1:1801,1], constant = FALSE, model = "cgarch",
                                         order = c(1,1), vreg = y[1:1801,2],
                                         distribution = "norm", init = "sample",
                                         sample_n = 100)
    local_mod_cgarch <- estimate(local_spec_cgarch)
    spec <- garch_modelspec(y[1:1000,1], constant = FALSE, model = "cgarch",
                            order = c(1,1), vreg = y[1:1000,2],
                            distribution = "norm", init = "sample",
                            sample_n = 100)
    spec$parmatrix <- copy(local_mod_cgarch$parmatrix)
    filt <- tsfilter(spec, y = y[1001:1801,1], newvreg = y[1001:1801,2])
    expect_equal(filt$sigma[1:1000], local_mod_cgarch$sigma[1:1000])
    expect_equal(filt$sigma[1:1801], local_mod_cgarch$sigma[1:1801])
    expect_equal(filt$permanent_component[1:1801], local_mod_cgarch$permanent_component[1:1801])
    expect_equal(filt$transitory_component[1:1801], local_mod_cgarch$transitory_component[1:1801])
    expect_equal(residuals(local_mod_cgarch, standardize = TRUE)[1:1801],
                 residuals(filt, standardize = T)[1:1801], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})

test_that("filter from spec (cgarch) append with initial conditions and multiplicative regressors: check results",{
    # see whether using same sample initialization with different data sizes
    # provides the same results
    local_spec_cgarch <- garch_modelspec(y[1:1801,1], constant = FALSE, model = "cgarch",
                                         order = c(1,1), vreg = y[1:1801,2], multiplicative = TRUE,
                                         distribution = "norm", init = "sample",
                                         sample_n = 100)
    local_mod_cgarch <- estimate(local_spec_cgarch)
    spec <- garch_modelspec(y[1:1000,1], constant = FALSE, model = "cgarch",
                            order = c(1,1), vreg = y[1:1000,2], multiplicative = TRUE,
                            distribution = "norm", init = "sample",
                            sample_n = 100)
    spec$parmatrix <- copy(local_mod_cgarch$parmatrix)
    filt <- tsfilter(spec, y = y[1001:1801,1], newvreg = y[1001:1801,2])
    expect_equal(filt$sigma[1:1000], local_mod_cgarch$sigma[1:1000])
    expect_equal(filt$sigma[1:1801], local_mod_cgarch$sigma[1:1801])
    expect_equal(filt$permanent_component[1:1801], local_mod_cgarch$permanent_component[1:1801])
    expect_equal(filt$transitory_component[1:1801], local_mod_cgarch$transitory_component[1:1801])
    expect_equal(residuals(local_mod_cgarch, standardize = TRUE)[1:1801],
                 residuals(filt, standardize = T)[1:1801], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})





test_that("filter from spec (egarch) append with initial conditions: check results",{
    # see whether using same sample initialization with different data sizes
    # provides the same results
    local_spec_egarch <- garch_modelspec(y[1:1801,1], constant = FALSE, model = "egarch",
                                         order = c(1,1), vreg = y[1:1801,2],
                                         distribution = "norm", init = "sample",
                                         sample_n = 100)
    local_mod_egarch <- estimate(local_spec_egarch)
    spec <- garch_modelspec(y[1:1000,1], constant = FALSE, model = "egarch",
                            order = c(1,1), vreg = y[1:1000,2],
                            distribution = "norm", init = "sample",
                            sample_n = 100)
    spec$parmatrix <- copy(local_mod_egarch$parmatrix)
    filt <- tsfilter(spec, y = y[1001:1801,1], newvreg = y[1001:1801,2])
    expect_equal(filt$sigma[1:1000], local_mod_egarch$sigma[1:1000])
    expect_equal(filt$sigma[1:1801], local_mod_egarch$sigma[1:1801])
    expect_equal(residuals(local_mod_egarch, standardize = TRUE)[1:1801],
                 residuals(filt, standardize = T)[1:1801], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})


test_that("filter from spec (gjrgarch) append with initial conditions: check results",{
    # see whether using same sample initialization with different data sizes
    # provides the same results
    local_spec_gjrgarch <- garch_modelspec(y[1:1801,1], constant = FALSE, model = "gjrgarch",
                                         order = c(1,1), vreg = y[1:1801,2],
                                         distribution = "norm", init = "sample",
                                         sample_n = 100)
    local_mod_gjrgarch <- estimate(local_spec_gjrgarch)
    spec <- garch_modelspec(y[1:1000,1], constant = FALSE, model = "gjrgarch",
                            order = c(1,1), vreg = y[1:1000,2],
                            distribution = "norm", init = "sample",
                            sample_n = 100)
    spec$parmatrix <- copy(local_mod_gjrgarch$parmatrix)
    filt <- tsfilter(spec, y = y[1001:1801,1], newvreg = y[1001:1801,2])
    expect_equal(filt$sigma[1:1000], local_mod_gjrgarch$sigma[1:1000])
    expect_equal(filt$sigma[1:1801], local_mod_gjrgarch$sigma[1:1801])
    expect_equal(residuals(local_mod_gjrgarch, standardize = TRUE)[1:1801],
                 residuals(filt, standardize = T)[1:1801], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})


test_that("filter from spec (fgarch) append with initial conditions: check results",{
    # see whether using same sample initialization with different data sizes
    # provides the same results
    local_spec_fgarch <- garch_modelspec(y[1:1801,1], constant = FALSE, model = "fgarch",
                                           order = c(1,1), vreg = y[1:1801,2],
                                           distribution = "std", init = "sample",
                                           sample_n = 100)
    local_mod_fgarch <- estimate(local_spec_fgarch)
    spec <- garch_modelspec(y[1:1000,1], constant = FALSE, model = "fgarch",
                            order = c(1,1), vreg = y[1:1000,2],
                            distribution = "std", init = "sample",
                            sample_n = 100)
    spec$parmatrix <- copy(local_mod_fgarch$parmatrix)
    filt <- tsfilter(spec, y = y[1001:1801,1], newvreg = y[1001:1801,2])
    expect_equal(filt$sigma[1:1000], local_mod_fgarch$sigma[1:1000])
    expect_equal(filt$sigma[1:1801], local_mod_fgarch$sigma[1:1801])
    expect_equal(residuals(local_mod_fgarch, standardize = TRUE)[1:1801],
                 residuals(filt, standardize = T)[1:1801], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})
