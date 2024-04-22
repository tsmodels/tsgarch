data(dmbp)
suppressMessages(suppressWarnings(library(xts)))
y <- xts(dmbp, as.Date(1:nrow(dmbp), origin = "1970-01-01"))

global_spec_garch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                               order = c(1,1), vreg = y[1:1800,2],
                               distribution = "norm")
global_mod_garch <- estimate(global_spec_garch)

global_spec_cgarch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "cgarch",
                                     order = c(1,1), vreg = y[1:1800,2],
                                     distribution = "norm")
global_mod_cgarch <- estimate(global_spec_cgarch)

global_spec_gjrgarch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "gjrgarch",
                                      order = c(1,1), vreg = y[1:1800,2],
                                      distribution = "norm")
global_mod_gjrgarch <- estimate(global_spec_gjrgarch)


global_spec_fgarch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "fgarch",
                                        order = c(1,1), vreg = y[1:1800,2],
                                        distribution = "norm")
global_mod_fgarch <- estimate(global_spec_fgarch)

global_spec_aparch <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "aparch",
                                      order = c(1,1), vreg = y[1:1800,2],
                                      distribution = "norm")
global_mod_aparch <- estimate(global_spec_aparch)


global_spec_garch_jsu <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                     order = c(1,1), vreg = y[1:1800,2],
                                     distribution = "jsu")
global_mod_garch_jsu <- estimate(global_spec_garch_jsu)


upper_trimmed_mean <- function(x, trim = 0) {
    x <- sort(x)
    mean(x[1:floor(length(x) * (1 - trim))])
}

