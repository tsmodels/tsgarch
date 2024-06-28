test_that("filter from estimate (garch) with tmb: check results",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                         order = c(1,1), vreg = y[1:1800,2],
                                         distribution = "norm")
    mod <- estimate(spec, keep_tmb = TRUE)
    lik_tmb <- mod$tmb$report(par = unname(coef(mod)))$ll_vector
    expect_equal(-log(lik_tmb), mod$lik_vector)
})

test_that("filter from estimate (gjrgarch) with tmb: check results",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch",
                            order = c(1,1), vreg = y[1:1800,2],
                            distribution = "norm")
    mod <- estimate(spec, keep_tmb = TRUE)
    lik_tmb <- mod$tmb$report(par = unname(coef(mod)))$ll_vector
    expect_equal(-log(lik_tmb), mod$lik_vector)
})
