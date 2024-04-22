test_that("gjr(1,1) newsimpact",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch",
                            order = c(1,1), distribution = "norm")
    mod <- estimate(spec)
    expect_no_error(newsimpact(mod))
})


