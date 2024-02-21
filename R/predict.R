# check multiplicative for other models
.predict_garch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = 1000, bootstrap = FALSE, forc_dates = NULL, init_states = NULL, ...)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))

    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                             new_regressors = newvreg, xi = object$parmatrix[group == "xi"]$value,
                             h = h, include_regressors = object$spec$vreg$include_vreg,
                             maxpq = maxpq,
                             regressor_argument = "newvreg")
    mu <- object$parmatrix[group == "mu"]$value
    omega <- omega(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    distribution <- object$parmatrix[group == "distribution"]$value
    xi = object$parmatrix[group == "xi"]$value
    constant <- omega + v
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    res <- c(tail(as.numeric(residuals(object)), maxpq), rep(0, h))
    sigma_sqr <- c(tail(as.numeric(sigma(object)^2), maxpq), rep(0, h))
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_sqr[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                sigma_sqr[i] <- sigma_sqr[i] + beta[j] * sigma_sqr[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- sigma_sqr[(i - j)]
                } else {
                    s <- res[(i - j)]^2
                }
                sigma_sqr[i] <- sigma_sqr[i] + alpha[j] * s
            }
        }
    }
    sigma <- sqrt(sigma_sqr)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    y <- y + mu
    if (bootstrap & h > 1) {
        out <- garch_boostrap(object, h, nsim, vreg = tail(v, h))
        sigma_sim <- out$sigma
        colnames(sigma_sim) <- as.character(forc_dates)
        class(sigma_sim) <- "tsmodel.distribution"
        attr(sigma_sim, "date_class") <- "Date"
        series_sim <- out$series
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    } else {
        series_sim <- do.call(cbind, lapply(1:h, function(i){
            rdist(distribution = object$spec$distribution, n = nsim, mu = mu, sigma = sigma[i], skew = distribution[1], shape = distribution[2], lambda = distribution[3])
        }))
        sigma_sim <- NULL
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    }
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = series_sim, sigma_sim  = sigma_sim, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}


.predict_egarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, bootstrap = FALSE, forc_dates = NULL, init_states = NULL, ...)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))

    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                             new_regressors = newvreg, xi = object$parmatrix[group == "xi"]$value,
                             h = h, include_regressors = object$spec$vreg$include_vreg,
                             maxpq = maxpq,
                             regressor_argument = "newvreg")
    mu <- object$parmatrix[group == "mu"]$value
    omega <- omega(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    distribution <- object$parmatrix[group == "distribution"]$value
    xi = object$parmatrix[group == "xi"]$value
    constant <- omega + v
    kappa <- object$kappa
    res <- c(tail(as.numeric(residuals(object)), maxpq), rep(0, h))
    std_res <- c(tail(as.numeric(residuals(object, standardize = TRUE)), maxpq), rep(0, h))
    sigma_sqr <- c(tail(as.numeric(sigma(object)^2), maxpq), rep(0, h))
    log_sigma_sqr <- c(tail(log(as.numeric(sigma(object)^2)), maxpq), rep(0, h))
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        log_sigma_sqr[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                log_sigma_sqr[i] <- log_sigma_sqr[i] + beta[j] * log_sigma_sqr[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- 0
                } else {
                    s <- alpha[j] * std_res[i - j] + gamma[j] * (abs(std_res[i - j]) - kappa)
                }
                log_sigma_sqr[i] <- log_sigma_sqr[i] + s
            }
        }
    }
    sigma <- sqrt(exp(log_sigma_sqr))
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    if (bootstrap & h > 1) {
        out <- garch_boostrap(object, h, nsim, vreg = tail(v, h))
        sigma_sim <- out$sigma
        colnames(sigma_sim) <- as.character(forc_dates)
        class(sigma_sim) <- "tsmodel.distribution"
        attr(sigma_sim, "date_class") <- "Date"
        series_sim <- out$series
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    } else {
        series_sim <- do.call(cbind, lapply(1:h, function(i){
            rdist(distribution = object$spec$distribution, n = nsim, mu = mu, sigma = sigma[i], skew = distribution[1], shape = distribution[2], lambda = distribution[3])
        }))
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
        sigma_sim <- NULL
    }
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = series_sim, sigma_sim  = sigma_sim, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}



.predict_aparch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, bootstrap = FALSE, forc_dates = NULL, init_states = NULL, ...)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))

    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                                        new_regressors = newvreg, xi = object$parmatrix[group == "xi"]$value,
                                        h = h, include_regressors = object$spec$vreg$include_vreg,
                                        maxpq = maxpq,
                                        regressor_argument = "newvreg")
    mu <- object$parmatrix[group == "mu"]$value
    omega <- omega(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    distribution <- object$parmatrix[group == "distribution"]$value
    xi = object$parmatrix[group == "xi"]$value
    constant <- omega + v
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    kappa <- object$kappa
    res <- c(tail(as.numeric(residuals(object)), maxpq), rep(0, h))
    sigma <- c(tail(as.numeric(sigma(object)), maxpq), rep(0, h))
    power_sigma <- c(tail(as.numeric(sigma(object)^delta), maxpq), rep(0, h))
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        power_sigma[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                power_sigma[i] <- power_sigma[i] + beta[j] * power_sigma[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- kappa[j] * power_sigma[i - j]
                } else {
                    s <- (abs(res[i - j]) - gamma[j] * res[i - j])^delta
                }
                power_sigma[i] <- power_sigma[i] + alpha[j] * s
            }
        }
    }
    sigma <- power_sigma^(1/delta)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    y <- y + mu
    if (bootstrap & h > 1) {
        out <- garch_boostrap(object, h, nsim, vreg = tail(v, h))
        sigma_sim <- out$sigma
        colnames(sigma_sim) <- as.character(forc_dates)
        class(sigma_sim) <- "tsmodel.distribution"
        attr(sigma_sim, "date_class") <- "Date"
        series_sim <- out$series
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    } else {
        series_sim <- do.call(cbind, lapply(1:h, function(i){
            rdist(distribution = object$spec$distribution, n = nsim, mu = mu, sigma = sigma[i], skew = distribution[1], shape = distribution[2], lambda = distribution[3])
        }))
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
        sigma_sim <- NULL
    }
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = series_sim, sigma_sim  = sigma_sim, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}

.predict_gjrgarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, bootstrap = FALSE, forc_dates = NULL, init_states = NULL, ...)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))

    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                                        new_regressors = newvreg, xi = object$parmatrix[group == "xi"]$value,
                                        h = h, include_regressors = object$spec$vreg$include_vreg,
                                        maxpq = maxpq,
                                        regressor_argument = "newvreg")
    mu <- object$parmatrix[group == "mu"]$value
    omega <- omega(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    distribution <- object$parmatrix[group == "distribution"]$value
    xi = object$parmatrix[group == "xi"]$value
    constant <- omega + v
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    kappa <- object$kappa
    res <- c(tail(as.numeric(residuals(object)), maxpq), rep(0, h))
    sigma <- c(tail(as.numeric(sigma(object)), maxpq), rep(0, h))
    sigma_squared <- c(tail(as.numeric(sigma(object)^2), maxpq), rep(0, h))
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_squared[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                sigma_squared[i] <- sigma_squared[i] + beta[j] * sigma_squared[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s1 <- sigma_squared[i - j]
                    s2 <- gamma[j] * kappa * sigma_squared[i - j]
                } else {
                    s1 <- res[i - j]^2
                    s2 <- gamma[j] * (res[i - j]^2 * (res[i - j] <= 0))
                }
                sigma_squared[i] <- sigma_squared[i] + alpha[j] * s1 + s2
            }
        }
    }
    sigma <- sqrt(sigma_squared)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]

    y <- y + mu
    if (bootstrap & h > 1) {
        out <- garch_boostrap(object, h, nsim, vreg = tail(v, h))
        sigma_sim <- out$sigma
        colnames(sigma_sim) <- as.character(forc_dates)
        class(sigma_sim) <- "tsmodel.distribution"
        attr(sigma_sim, "date_class") <- "Date"
        series_sim <- out$series
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    } else {
        series_sim <- do.call(cbind, lapply(1:h, function(i){
            rdist(distribution = object$spec$distribution, n = nsim, mu = mu, sigma = sigma[i], skew = distribution[1], shape = distribution[2], lambda = distribution[3])
        }))
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
        sigma_sim <- NULL
    }
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = series_sim, sigma_sim  = sigma_sim, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}

.predict_fgarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, bootstrap = FALSE, forc_dates = NULL, init_states = NULL, ...)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))

    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                                        new_regressors = newvreg, xi = object$parmatrix[group == "xi"]$value,
                                        h = h, include_regressors = object$spec$vreg$include_vreg,
                                        maxpq = maxpq,
                                        regressor_argument = "newvreg")
    mu <- object$parmatrix[group == "mu"]$value
    omega <- omega(object)
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    eta <- object$parmatrix[group == "eta"]$value
    beta <- object$parmatrix[group == "beta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    distribution <- object$parmatrix[group == "distribution"]$value
    xi = object$parmatrix[group == "xi"]$value
    constant <- omega + v
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    kappa <- object$kappa
    res <- c(tail(as.numeric(residuals(object)), maxpq), rep(0, h))
    sigma <- c(tail(as.numeric(sigma(object)), maxpq), rep(0, h))
    std_res <- c(tail(as.numeric(residuals(object, standardize = TRUE)), maxpq), rep(0, h))
    power_sigma <- c(tail(as.numeric(sigma(object)^delta), maxpq), rep(0, h))
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        power_sigma[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                power_sigma[i] <- power_sigma[i] + beta[j] * power_sigma[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- power_sigma[i - j] * kappa[j]
                } else {
                    s <-  power_sigma[i - j] * (abs(std_res[i - j] - eta[j]) - gamma[j] * (std_res[i - j] - eta[j]))^delta
                }
                power_sigma[i] <- power_sigma[i] + alpha[j] * s
            }
        }
    }
    sigma <- power_sigma^(1/delta)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]

    y <- y + mu
    if (bootstrap & h > 1) {
        out <- garch_boostrap(object, h, nsim, vreg = tail(v, h))
        sigma_sim <- out$sigma
        colnames(sigma_sim) <- as.character(forc_dates)
        class(sigma_sim) <- "tsmodel.distribution"
        attr(sigma_sim, "date_class") <- "Date"
        series_sim <- out$series
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    } else {
        series_sim <- do.call(cbind, lapply(1:h, function(i){
            rdist(distribution = object$spec$distribution, n = nsim, mu = mu, sigma = sigma[i], skew = distribution[1], shape = distribution[2], lambda = distribution[3])
        }))
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
        sigma_sim <- NULL
    }
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = series_sim, sigma_sim  = sigma_sim, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}

.predict_cgarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, bootstrap = FALSE, forc_dates = NULL, init_states = NULL, ...)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))

    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                                        new_regressors = newvreg, xi = object$parmatrix[group == "xi"]$value,
                                        h = h, include_regressors = object$spec$vreg$include_vreg,
                                        maxpq = maxpq,
                                        regressor_argument = "newvreg")
    mu <- object$parmatrix[group == "mu"]$value
    omega <- omega(object)
    rho <- object$parmatrix[group == "rho"]$value
    phi <- object$parmatrix[group == "phi"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    distribution <- object$parmatrix[group == "distribution"]$value
    xi = object$parmatrix[group == "xi"]$value
    constant <- omega + v
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    res <- c(tail(as.numeric(residuals(object)), maxpq), rep(0, h))
    sigma_sqr <- c(tail(as.numeric(sigma(object)^2), maxpq), rep(0, h))
    permanent_component <- c(tail(as.numeric(object$permanent_component), maxpq), rep(0, h))
    y <- rep(0, h)
    k <- 1
    for (i in (maxpq + 1):(h + maxpq)) {
        if (k == 1) {
            permanent_component[i] <- constant[i] + rho * permanent_component[i - 1] + phi * (res[i - 1]^2 - sigma_sqr[i - 1])
            sigma_sqr[i] <- permanent_component[i]
            if (garch_order[2] > 0) {
                for (j in 1:garch_order[2]) {
                    sigma_sqr[i] <- sigma_sqr[i] + beta[j] * (sigma_sqr[i - j] - permanent_component[i - j])
                }
            }
            if (garch_order[1] > 0) {
                for (j in 1:garch_order[1]) {
                    if ((i - maxpq) > j) {
                        s <- sigma_sqr[(i - j)]
                    } else {
                        s <- res[(i - j)]^2
                    }
                    sigma_sqr[i] <- sigma_sqr[i] + alpha[j] * (s - permanent_component[i - j])
                }
            }
        } else {
            permanent_component[i] <- constant[i]/(1 - rho) + rho^k * (permanent_component[maxpq] - constant[i]/(1 - rho))
            p <- sum(beta) + sum(alpha)
            sigma_sqr[i] <- permanent_component[i]  +  p^k * (sigma_sqr[maxpq] - permanent_component[maxpq])
        }
        k <- k + 1
    }
    sigma <- sqrt(sigma_sqr)
    if (maxpq > 0) {
        sigma <- sigma[-seq_len(maxpq)]
        permanent_component <- permanent_component[-seq_len(maxpq)]
    }

    y <- y + mu
    if (bootstrap & h > 1) {
        out <- garch_boostrap(object, h, nsim, vreg = tail(v, h))
        sigma_sim <- out$sigma
        colnames(sigma_sim) <- as.character(forc_dates)
        class(sigma_sim) <- "tsmodel.distribution"
        attr(sigma_sim, "date_class") <- "Date"
        series_sim <- out$series
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
    } else {
        series_sim <- do.call(cbind, lapply(1:h, function(i){
            rdist(distribution = object$spec$distribution, n = nsim, mu = mu, sigma = sigma[i], skew = distribution[1], shape = distribution[2], lambda = distribution[3])
        }))
        colnames(series_sim) <- as.character(forc_dates)
        class(series_sim) <- "tsmodel.distribution"
        attr(series_sim, "date_class") <- "Date"
        sigma_sim <- NULL
    }
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    permanent_component <- xts(permanent_component, forc_dates)
    transitory_component <- sigma - permanent_component
    L <- list(original_series = object$spec$target$y, distribution = series_sim,
              sigma_sim  = sigma_sim, spec = object$spec, sigma = sigma,
              permanent_component = permanent_component,
              transitory_component = transitory_component,
              mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}


garch_boostrap <- function(object, h, nsim, vreg)
{
    z <- as.numeric(residuals(object, standardize = TRUE))
    zsim <- do.call(rbind, lapply(1:nsim, function(i){
        sample(z, h, replace = TRUE)
    }))
    spec <- object$spec
    spec$parmatrix <- object$parmatrix
    maxpq <- max(spec$model$order)
    init_v <- tail(as.numeric(object$sigma), maxpq)^2
    init_z <- tail(z, maxpq)
    # if model == cgarch must provide a matrix for var_init
    if (object$spec$model$model == "cgarch") {
        init_v <- cbind(init_v, tail(object$permanent_component, maxpq))
    }
    b <- simulate(spec, h = h, nsim = nsim, var_init = init_v, innov = zsim, innov_init = init_z, vreg = vreg)
    return(b)
}
