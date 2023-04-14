initialize_parameters <- function(model = "garch", y, constant = 0.0,
                                  order = c(1,1), variance_targeting = FALSE,
                                  vreg = NULL, multiplicative = TRUE,
                                  init = c("unconditional","sample","backcast"),
                                  backcast_lambda = 0.7, sample_n = 10,
                                  distribution = "norm", ...)
{
    switch(model,
           "garch" = .parameters_garch(y = y, constant = constant,
                                       order = order,
                                       variance_targeting = variance_targeting,
                                       vreg = vreg,
                                       multiplicative = multiplicative,
                                       init = init,
                                       backcast_lambda = backcast_lambda,
                                       sample_n = sample_n,
                                       distribution = distribution),
           "egarch" = .parameters_egarch(y = y, constant = constant,
                                         order = order,
                                         variance_targeting = variance_targeting,
                                         vreg = vreg,
                                         multiplicative = multiplicative,
                                         init = init,
                                         backcast_lambda = backcast_lambda,
                                         sample_n = sample_n,
                                         distribution = distribution),
           "aparch" = .parameters_aparch(y = y, constant = constant,
                                         order = order,
                                         variance_targeting = variance_targeting,
                                         vreg = vreg,
                                         multiplicative = multiplicative,
                                         init = init,
                                         backcast_lambda = backcast_lambda,
                                         sample_n = sample_n,
                                         distribution = distribution),
           "gjrgarch" = .parameters_gjrgarch(y = y, constant = constant,
                                         order = order,
                                         variance_targeting = variance_targeting,
                                         vreg = vreg,
                                         multiplicative = multiplicative,
                                         init = init,
                                         backcast_lambda = backcast_lambda,
                                         sample_n = sample_n,
                                         distribution = distribution),
           "fgarch" = .parameters_fgarch(y = y, constant = constant,
                                             order = order,
                                             variance_targeting = variance_targeting,
                                             vreg = vreg,
                                             multiplicative = multiplicative,
                                             init = init,
                                             backcast_lambda = backcast_lambda,
                                             sample_n = sample_n,
                                             distribution = distribution),
           "cgarch" = .parameters_cgarch(y = y, constant = constant,
                                         order = order,
                                         variance_targeting = variance_targeting,
                                         vreg = vreg,
                                         multiplicative = multiplicative,
                                         init = init,
                                         backcast_lambda = backcast_lambda,
                                         sample_n = sample_n,
                                         distribution = distribution),
           "igarch" = .parameters_igarch(y = y, constant = constant,
                                         order = order,
                                         variance_targeting = variance_targeting,
                                         vreg = vreg,
                                         multiplicative = multiplicative,
                                         init = init,
                                         backcast_lambda = backcast_lambda,
                                         sample_n = sample_n,
                                         distribution = distribution),
           "ewma" = .parameters_ewma(y = y, constant = constant,
                                         order = order,
                                         variance_targeting = variance_targeting,
                                         vreg = vreg,
                                         multiplicative = multiplicative,
                                         init = init,
                                         backcast_lambda = backcast_lambda,
                                         sample_n = sample_n,
                                         distribution = distribution))
}

.parameters_garch <- function(y, constant = FALSE, order = c(1,1),
                              variance_targeting = FALSE,
                              vreg = NULL, multiplicative = TRUE,
                              init = c("unconditional","sample","backcast"),
                              backcast_lambda = 0.7, sample_n = 10,
                              distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100,
                            upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0),
                            scale = 1, group = "mu", equation = "[M]",
                            symbol = "\\mu")

    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    } else {
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2)) * 0.01]
    }
    if (!is.null(vreg)) {
        if (multiplicative) {
            parmatrix[parameter == "omega", lower := log(1e-12)]
            parmatrix[parameter == "omega", upper :=  log(upper)]
            parmatrix[parameter == "omega", value :=  log(value)]
        }
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.05/order[1], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "alpha",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0, lower = 0,
                                      upper = 1, estimate = 0, scale = 1,
                                      group = "beta", equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.9/order[2], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = 0, upper = 100, estimate = 1, scale = 1,
                              group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
        if (multiplicative) {
            vmatrix$lower <- -1000
            vmatrix$upper <- 1000
        }
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0, lower = 0,
                              upper = 1, estimate = 0, scale = 1, group = "xi",
                              equation = "[V]", symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}


.parameters_egarch <- function(y, constant = FALSE, order = c(1,1),
                               variance_targeting = FALSE,
                               vreg = NULL, multiplicative = TRUE,
                               init = c("unconditional","sample","backcast"),
                               backcast_lambda = 0.7, sample_n = 10,
                               distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100, upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0), scale = 1,
                            group = "mu", equation = "[M]", symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = log(var_y) * 0.01,
                                  lower = -10, upper = 10, estimate = 1,
                                  scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = -10, upper = 10, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "gamma1", value = 0,
                                      lower = -10, upper = 10, estimate = 0,
                                      scale = 1, group = "gamma",
                                      equation = "[ARCH]", symbol = "\\gamma_1"))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.1/order[1], lower = -10,
                                      upper = 10, estimate = 1, scale = 1,
                                      group = "alpha", equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
        parmatrix <- rbind(parmatrix,
              data.table("parameter" = paste0("gamma",1:order[1]),
                         value = 0.1/order[1], lower = -10, upper = 10,
                         estimate = 1, scale = 1, group = "gamma",
                         equation = "[ARCH]",
                         symbol = paste0("\\gamma_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0,
                                      lower = -1 + 1e-12, upper = 1 - 1e-12,
                                      estimate = 0, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.9/order[2], lower = -1 + 1e-12,
                                      upper =  1 - 1e-12, estimate = 1,
                                      scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = -1000, upper = 1000, estimate = 1,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0,
                              lower = -1000, upper = 1000, estimate = 0,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}

.parameters_aparch <- function(y, constant = FALSE, order = c(1,1),
                               variance_targeting = FALSE,
                               vreg = NULL, multiplicative = TRUE,
                               init = c("unconditional","sample","backcast"),
                               backcast_lambda = 0.7, sample_n = 10,
                               distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    delta <- 2
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100, upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0), scale = 1,
                            group = "mu", equation = "[M]", symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1 - 1e-3, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "gamma1", value = 0,
                                      lower = -1 + 1e-3, upper = 1 - 1e-3, estimate = 0,
                                      scale = 1, group = "gamma",
                                      equation = "[ARCH]", symbol = "\\gamma_1"))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.1/order[1], lower = 0,
                                      upper = 1 - 1e-3, estimate = 1, scale = 1,
                                      group = "alpha", equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("gamma",1:order[1]),
                                      value = 0.05/order[1], lower = -1 + 1e-3,
                                      upper = 1 - 1e-3,
                                      estimate = 1, scale = 1, group = "gamma",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\gamma_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0,
                                      lower = 0, upper = 1 - 1e-3,
                                      estimate = 0, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.8/order[2], lower = 0,
                                      upper =  1 - 1e-3, estimate = 1,
                                      scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    parmatrix <- rbind(parmatrix, data.table("parameter" = "delta", value = delta,
                                             lower = 0.01, upper = 3.5,
                                             estimate = 1,
                                             scale = 1, group = "delta",
                                             equation = "[GARCH]",
                                             symbol = "\\delta"))
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = -1000, upper = 1000, estimate = 1,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0,
                              lower = -1000, upper = 1000, estimate = 0,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}


.parameters_gjrgarch <- function(y, constant = FALSE, order = c(1,1),
                               variance_targeting = FALSE,
                               vreg = NULL, multiplicative = TRUE,
                               init = c("unconditional","sample","backcast"),
                               backcast_lambda = 0.7, sample_n = 10,
                               distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100, upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0), scale = 1,
                            group = "mu", equation = "[M]", symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1 - 1e-12, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "gamma1", value = 0,
                                      lower = -1 + 1e-12, upper = 1 - 1e-12, estimate = 0,
                                      scale = 1, group = "gamma",
                                      equation = "[ARCH]", symbol = "\\gamma_1"))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.1/order[1], lower = 0,
                                      upper = 1 - 1e-12, estimate = 1, scale = 1,
                                      group = "alpha", equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("gamma",1:order[1]),
                                      value = 0.05/order[1], lower = -1 + 1e-12,
                                      upper = 1 - 1e-12,
                                      estimate = 1, scale = 1, group = "gamma",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\gamma_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0,
                                      lower = 0, upper = 1 - 1e-12,
                                      estimate = 0, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.8/order[2], lower = 0,
                                      upper =  1 - 1e-12, estimate = 1,
                                      scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = -1000, upper = 1000, estimate = 1,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0,
                              lower = -1000, upper = 1000, estimate = 0,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}

.parameters_fgarch <- function(y, constant = FALSE, order = c(1,1),
                               variance_targeting = FALSE,
                               vreg = NULL, multiplicative = TRUE,
                               init = c("unconditional","sample","backcast"),
                               backcast_lambda = 0.7, sample_n = 10,
                               distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    delta <- 2
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100, upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0), scale = 1,
                            group = "mu", equation = "[M]", symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1 - 1e-12, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "gamma1", value = 0,
                                      lower = -1 + 1e-12, upper = 1 - 1e-12, estimate = 0,
                                      scale = 1, group = "gamma",
                                      equation = "[ARCH]", symbol = "\\gamma_1"))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "eta1", value = 0,
                                      lower = -10, upper = 10, estimate = 0,
                                      scale = 1, group = "eta",
                                      equation = "[ARCH]", symbol = "\\eta_1"))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.1/order[1], lower = 0,
                                      upper = 1 - 1e-12, estimate = 1, scale = 1,
                                      group = "alpha", equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("gamma",1:order[1]),
                                      value = 0.05/order[1], lower = -1 + 1e-12,
                                      upper = 1 - 1e-12,
                                      estimate = 1, scale = 1, group = "gamma",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\gamma_",1:order[1])))
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("eta",1:order[1]),
                                      value = 0.05/order[1], lower = -10,
                                      upper = 10,
                                      estimate = 1, scale = 1, group = "eta",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\eta_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0,
                                      lower = 0, upper = 1 - 1e-12,
                                      estimate = 0, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.8/order[2], lower = 0,
                                      upper =  1 - 1e-12, estimate = 1,
                                      scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    parmatrix <- rbind(parmatrix, data.table("parameter" = "delta", value = delta,
                                             lower = 0.01, upper = 3.5,
                                             estimate = 1,
                                             scale = 1, group = "delta",
                                             equation = "[GARCH]",
                                             symbol = "\\delta"))
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = -1000, upper = 1000, estimate = 1,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0,
                              lower = -1000, upper = 1000, estimate = 0,
                              scale = 1, group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}

.parameters_cgarch <- function(y, constant = FALSE, order = c(1,1),
                              variance_targeting = FALSE,
                              vreg = NULL, multiplicative = TRUE,
                              init = c("unconditional","sample","backcast"),
                              backcast_lambda = 0.7, sample_n = 10,
                              distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100,
                            upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0),
                            scale = 1, group = "mu", equation = "[M]",
                            symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    } else {
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2)) * 0.01]
    }
    if (!is.null(vreg)) {
        if (multiplicative) {
            parmatrix[parameter == "omega", lower := log(1e-12)]
            parmatrix[parameter == "omega", upper :=  log(upper)]
            parmatrix[parameter == "omega", value :=  log(value)]
        }
    }
    # permanent component parameters \\rho and \\phi
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "rho", value = 0.98,
                                  lower = 0, upper = 1 - 1e-12,
                                  estimate = 1, scale = 1, group = "rho",
                                  equation = "[V]", symbol = "\\rho"))

    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "phi", value = 0.025/order[1],
                                  lower = 0 + 1e-12, upper = 1 - 1e-12,
                                  estimate = 1, scale = 1, group = "phi",
                                  equation = "[V]", symbol = "\\phi"))

    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.05/order[1], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "alpha",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0, lower = 0,
                                      upper = 1, estimate = 0, scale = 1,
                                      group = "beta", equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.9/order[2], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = 0, upper = 100, estimate = 1, scale = 1,
                              group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
        if (multiplicative) {
            vmatrix$lower <- -1000
            vmatrix$upper <- 1000
        }
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0, lower = 0,
                              upper = 1, estimate = 0, scale = 1, group = "xi",
                              equation = "[V]", symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}


.parameters_igarch <- function(y, constant = FALSE, order = c(1,1),
                              variance_targeting = FALSE,
                              vreg = NULL, multiplicative = TRUE,
                              init = c("unconditional","sample","backcast"),
                              backcast_lambda = 0.7, sample_n = 10,
                              distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100,
                            upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0),
                            scale = 1, group = "mu", equation = "[M]",
                            symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    if (variance_targeting) {
        parmatrix[parameter == "omega", estimate := 0]
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2))]
    } else {
        parmatrix[parameter == "omega", value := as.numeric(mean((y - mu)^2)) * 0.01]
    }
    if (!is.null(vreg)) {
        if (multiplicative) {
            parmatrix[parameter == "omega", lower := log(1e-12)]
            parmatrix[parameter == "omega", upper :=  log(upper)]
            parmatrix[parameter == "omega", value :=  log(value)]
        }
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.05/order[1], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "alpha",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0, lower = 0,
                                      upper = 1, estimate = 0, scale = 1,
                                      group = "beta", equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.9/order[2], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    # now check for summation to 1 and other conditions
    if (order[2] == 0 &
        length(parmatrix[group == "alpha" & estimate == 0]$value) == order[1] &
        sum(parmatrix[group == "alpha" & estimate == 0]$value) < 1) {
        stop("\ncannot estimate an igarch model with GARCH(p=0) and GARCH(q) fixed but
                 not summing to 1")
    }
    if (order[1] == 0 &
        length(parmatrix[group == "beta" & estimate == 0]$value) == order[2] &
        sum(parmatrix[group == "beta" & estimate == 0]$value) < 1) {
        stop("\ncannot estimate an igarch model with GARCH(q=0) and GARCH(p) fixed but
                 not summing to 1")
    }
    sum_alpha <- sum(parmatrix[group == "alpha"]$value)
    parmatrix[group == "beta", value := (1 - sum_alpha)/order[2]]
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = 0, upper = 100, estimate = 1, scale = 1,
                              group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
        if (multiplicative) {
            vmatrix$lower <- -1000
            vmatrix$upper <- 1000
        }
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0, lower = 0,
                              upper = 1, estimate = 0, scale = 1, group = "xi",
                              equation = "[V]", symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}

.parameters_ewma <- function(y, constant = FALSE, order = c(1,1),
                               variance_targeting = FALSE,
                               vreg = NULL, multiplicative = TRUE,
                               init = c("unconditional","sample","backcast"),
                               backcast_lambda = 0.7, sample_n = 10,
                               distribution = "norm", ...)
{
    parameter <- value <- group <- lower <- upper <- NULL
    y <- as.numeric(y)
    maxpq <- max(order)
    n <- NROW(y)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    var_y <- var(y, na.rm = T)
    parmatrix <- data.table("parameter" = "mu", value = mu,
                            lower =  -1.0 * abs(mu) * 100,
                            upper = abs(mu) * 100,
                            estimate = ifelse(constant, 1, 0),
                            scale = 1, group = "mu", equation = "[M]",
                            symbol = "\\mu")
    parmatrix <- rbind(parmatrix,
                       data.table("parameter" = "omega", value = var_y * 0.01,
                                  lower = 1e-12, upper = var_y/0.01,
                                  estimate = 1, scale = 1, group = "omega",
                                  equation = "[V]", symbol = "\\omega"))
    parmatrix[parameter == "omega", estimate := 0]
    parmatrix[parameter == "omega", value := 0]
    if (!is.null(vreg)) {
        if (multiplicative) {
            parmatrix[parameter == "omega", lower := log(1e-12)]
            parmatrix[parameter == "omega", upper :=  log(upper)]
            parmatrix[parameter == "omega", value :=  log(value)]
        }
    }
    if (order[1] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "alpha1", value = 0,
                                      lower = 0, upper = 1, estimate = 0,
                                      scale = 1, group = "alpha",
                                      equation = "[ARCH]", symbol = "\\alpha_1"))
    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("alpha",1:order[1]),
                                      value = 0.05/order[1], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "alpha",
                                      equation = "[ARCH]",
                                      symbol = paste0("\\alpha_",1:order[1])))
    }
    if (order[2] == 0) {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = "beta", value = 0, lower = 0,
                                      upper = 1, estimate = 0, scale = 1,
                                      group = "beta", equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1)))

    } else {
        parmatrix <- rbind(parmatrix,
                           data.table("parameter" = paste0("beta",1:order[2]),
                                      value = 0.9/order[2], lower = 0, upper = 1,
                                      estimate = 1, scale = 1, group = "beta",
                                      equation = "[GARCH]",
                                      symbol = paste0("\\beta_",1:order[2])))

    }
    # now check for summation to 1 and other conditions
    if (order[2] == 0 &
        length(parmatrix[group == "alpha" & estimate == 0]$value) == order[1] &
        sum(parmatrix[group == "alpha" & estimate == 0]$value) < 1) {
        stop("\ncannot estimate an igarch model with GARCH(p=0) and GARCH(q) fixed but
                 not summing to 1")
    }
    if (order[1] == 0 &
        length(parmatrix[group == "beta" & estimate == 0]$value) == order[2] &
        sum(parmatrix[group == "beta" & estimate == 0]$value) < 1) {
        stop("\ncannot estimate an igarch model with GARCH(q=0) and GARCH(p) fixed but
                 not summing to 1")
    }
    sum_alpha <- sum(parmatrix[group == "alpha"]$value)
    parmatrix[group == "beta", value := (1 - sum_alpha)/order[2]]
    if (!is.null(vreg)) {
        m <- NCOL(vreg)
        vmatrix <- data.table("parameter" = paste0("xi",1:m), value = 1,
                              lower = 0, upper = 100, estimate = 1, scale = 1,
                              group = "xi", equation = "[V]",
                              symbol = paste0("\\xi_",1:m))
        if (multiplicative) {
            vmatrix$lower <- -1000
            vmatrix$upper <- 1000
        }
    } else {
        vmatrix <- data.table("parameter" = paste0("xi",1), value = 0, lower = 0,
                              upper = 1, estimate = 0, scale = 1, group = "xi",
                              equation = "[V]", symbol = paste0("\\xi_",1))
    }
    parmatrix <- rbind(parmatrix, vmatrix)
    dmatrix <- distribution_parameters(distribution)
    parmatrix <- rbind(parmatrix, dmatrix)
    parmatrix[,estimate := as.integer(estimate)]
    return(parmatrix)
}
