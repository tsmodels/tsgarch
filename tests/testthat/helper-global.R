data(dmbp)
suppressMessages(suppressWarnings(library(xts)))
y <- xts(dmbp, as.Date(1:nrow(dmbp), origin = "1970-01-01"))
global_spec_norm <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                               order = c(1,1), vreg = y[1:1800,2],
                               distribution = "norm")
global_mod_norm <- estimate(global_spec_norm)

global_spec_ghst <- garch_modelspec(y[1:1800,1], constant = TRUE, model = "garch",
                                    order = c(1,1), vreg = y[1:1800,2], multiplicative = TRUE,
                                    distribution = "ghst")
global_spec_ghst$parmatrix[parameter == "shape", lower := 4.1]
global_mod_ghst <- estimate(global_spec_ghst)
