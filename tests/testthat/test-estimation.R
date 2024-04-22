test_that("estimate: garch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1076.574, tolerance = 1e-3)
})

test_that("estimate: garch(1,1)-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.136169]
    spec$parmatrix[parameter == "beta1", value := 0.823269]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1076.574, tolerance = 1e-3)
})

test_that("estimate: egarch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch", order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1073.671, tolerance = 1e-3)
})


test_that("estimate: egarch(1,1)-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "egarch", order = c(1,1))
    spec$parmatrix[parameter == "beta1", value := 0.9174076153]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1073.671, tolerance = 1e-3)
})

test_that("estimate: gjrgarch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch", order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1076.041, tolerance = 1e-3)
})

test_that("estimate: gjrgarch(1,1)-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.123969]
    spec$parmatrix[parameter == "beta1", value := 0.816408]
    spec$parmatrix[parameter == "gamma1", value := 0.028495]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    spec$parmatrix[parameter == "gamma1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1076.041, tolerance = 1e-3)
})

test_that("estimate: aparch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch", order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1074.1935, tolerance = 1e-3)
})

test_that("estimate: aparch(1,1)-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.157982]
    spec$parmatrix[parameter == "beta1", value := 0.807856]
    spec$parmatrix[parameter == "gamma1", value := 0.100214]
    spec$parmatrix[parameter == "delta", value := 1.439003]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    spec$parmatrix[parameter == "gamma1", estimate := 0]
    spec$parmatrix[parameter == "delta", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1074.1935, tolerance = 1e-3)
})

test_that("estimate: fgarch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch", order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1073.5190, tolerance = 1e-3)
})

test_that("estimate: fgarch(1,1)-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.154864]
    spec$parmatrix[parameter == "beta1", value := 0.810884]
    spec$parmatrix[parameter == "gamma1", value := -0.055386]
    spec$parmatrix[parameter == "eta1", value := 0.230461]
    spec$parmatrix[parameter == "delta", value := 1.582117]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    spec$parmatrix[parameter == "gamma1", estimate := 0]
    spec$parmatrix[parameter == "delta", estimate := 0]
    spec$parmatrix[parameter == "eta1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1073.5190, tolerance = 1e-3)
})

test_that("estimate: fgarch(1,1)-jsu-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch", order = c(1,1), distribution = "jsu")
    spec$parmatrix[parameter == "alpha1", value := 0.154864]
    spec$parmatrix[parameter == "beta1", value := 0.810884]
    spec$parmatrix[parameter == "gamma1", value := -0.055386]
    spec$parmatrix[parameter == "eta1", value := 0.230461]
    spec$parmatrix[parameter == "delta", value := 1.582117]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    spec$parmatrix[parameter == "gamma1", estimate := 0]
    spec$parmatrix[parameter == "delta", estimate := 0]
    spec$parmatrix[parameter == "eta1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -978.4525, tolerance = 1e-3)
})


test_that("estimate: cgarch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch", order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1060.7082, tolerance = 1e-3)
})


test_that("estimate: cgarch(1,1)-partially fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.1528533]
    spec$parmatrix[parameter == "beta1", value := 0.5079146]
    spec$parmatrix[parameter == "rho", value := 0.9938057]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    spec$parmatrix[parameter == "rho", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1060.7082, tolerance = 1e-3)
})

test_that("estimate: cgarch(1,1)-fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.1528533]
    spec$parmatrix[parameter == "beta1", value := 0.5079146]
    spec$parmatrix[parameter == "rho", value := 0.9938057]
    spec$parmatrix[parameter == "phi", value := 0.0440527]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    spec$parmatrix[parameter == "rho", estimate := 0]
    spec$parmatrix[parameter == "phi", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1060.7082, tolerance = 1e-3)
})

test_that("estimate: igarch(1,1)",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "igarch", order = c(1,1))
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1082.6435, tolerance = 1e-3)
})

test_that("estimate: igarch(1,1)-partially fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "igarch", order = c(1,1))
    spec$parmatrix[parameter == "alpha1", value := 0.140900]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", value := 0.859100]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(as.numeric(logLik(mod)), -1082.6435, tolerance = 1e-3)
})

test_that("estimate: igarch(2,1)-partially fixed parameters",{
    spec <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "igarch", order = c(2,1))
    spec$parmatrix[parameter == "alpha1", value := 0.100900]
    spec$parmatrix[parameter == "alpha1", estimate := 0]
    spec$parmatrix[parameter == "beta1", value := 0.859100]
    spec$parmatrix[parameter == "beta1", estimate := 0]
    mod <- estimate(spec)
    expect_equal(persistence(mod), 1.0)
    expect_equal(as.numeric(logLik(mod)), -1088.12, tolerance = 1e-3)
})
