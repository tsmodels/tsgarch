# check multiplicative for other models
get_group_parameters <- function(object)
{
    group <- NULL
    model <- object$spec$model$model
    parlist <- list()
    parlist$mu <- object$parmatrix[group == "mu"]$value
    parlist$omega <- omega(object)
    parlist$alpha <- object$parmatrix[group == "alpha"]$value
    parlist$beta <- object$parmatrix[group == "beta"]$value
    parlist$xi <- object$parmatrix[group == "xi"]$value
    parlist$distribution <- object$parmatrix[group == "distribution"]$value
    if (model %in% c("egarch","gjrgarch")) {
        parlist$gamma <- object$parmatrix[group == "gamma"]$value
    } else if (model == "aparch") {
        parlist$gamma <- object$parmatrix[group == "gamma"]$value
        parlist$delta <- object$parmatrix[group == "delta"]$value
    } else if (model == "fgarch") {
        parlist$gamma <- object$parmatrix[group == "gamma"]$value
        parlist$eta <- object$parmatrix[group == "eta"]$value
        parlist$delta <- object$parmatrix[group == "delta"]$value
    } else if (model == "cgarch") {
        parlist$phi <- object$parmatrix[group == "phi"]$value
        parlist$rho <- object$parmatrix[group == "rho"]$value
    }
    return(parlist)
}

setup_prediction <- function(object, h, newxreg = NULL, newvreg = NULL, forc_dates = NULL)
{
    parameter <- group <- NULL
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    forc_dates <- .forecast_dates(forc_dates = forc_dates, h = h, sampling = object$spec$target$sampling,
                                  last_index = tail(object$spec$target$index, 1))
    model_parameters <- get_group_parameters(object)
    v <- .process_prediction_regressors(old_regressors = object$spec$vreg$vreg,
                                        new_regressors = newvreg, xi = model_parameters$xi,
                                        h = h, include_regressors = object$spec$vreg$include_vreg,
                                        maxpq = maxpq, regressor_argument = "newvreg")

    constant <- model_parameters$omega + v
    return(list(variance_regressors = v, model_parameters = model_parameters, forc_dates = forc_dates, garch_order = garch_order, maxpq = maxpq))
}


# arch_init
# GARCH/APARCH -> squared_residuals
# init = list(variance = NULL, residuals = NULL, permanent_component = NULL, transitory_component = NULL)
initialize_states <- function(object, init_states = list())
{
    model <- object$spec$model$model
    maxpq <- max(object$spec$model$order)
    # switch to [[]] to avoid partial matching (silly behavior in R)
    if (is.null(init_states[["variance"]])) {
        init_variance <- tail(as.numeric(sigma(object)^2), maxpq)
    } else {
        if (length(init_states[["variance"]]) != maxpq) stop(paste0("\ninit_states$variance must be of length max (p,q): ", maxpq))
        init_variance <- init_states[["variance"]]
    }
    if (is.null(init_states[["residuals"]])) {
        init_residuals <- tail(as.numeric(residuals(object)), maxpq)
    } else {
        if (length(init_states[["residuals"]]) != maxpq) stop(paste0("\ninit_states$residuals must be of length max (p,q): ", maxpq))
        init_residuals <- init_states[["residuals"]]
    }

    if (is.null(init_states[["std_residuals"]])) {
        init_std_residuals <- tail(as.numeric(residuals(object, standardize = TRUE)), maxpq)
    } else {
        if (length(init_states[["std_residuals"]]) != maxpq) stop(paste0("\ninit_states$std_residuals must be of length max (p,q): ", maxpq))
        init_std_residuals <- init_states[["std_residuals"]]
    }

    if (model == "cgarch") {
        if (is.null(init_states$permanent_component)) {
            init_permanent_component <- tail(object$permanent_component, maxpq)
        } else {
            if (length(init_states$permanent_component) != maxpq) stop(paste0("\ninit_states$permanent_component must be of length max (p,q): ", maxpq))
            init_permanent_component <- init_states$permanent_component
        }
        if (is.null(init_states$transitory_component)) {
            init_transitory_component <- tail(object$transitory_component, maxpq)
        } else {
            if (length(init_states$transitory_component) != maxpq) stop(paste0("\ninit_states$transitory_component must be of length max (p,q): ", maxpq))
            init_transitory_component <- init_states$transitory_component
        }
    } else {
        init_permanent_component <- NULL
        init_transitory_component <- NULL
    }
    return(list(variance = init_variance, residuals = init_residuals, std_residuals = init_std_residuals, permanent_component = init_permanent_component,
                transitory_component = init_transitory_component))
}


simulated_distribution <- function(object, sigma, h = 1, nsim = 1,
                                   sim_method = c("parametric","bootstrap"),
                                   block = 1, model_parameters,
                                   vreg = NULL, forc_dates = NULL,
                                   bootstrap = FALSE, seed = NULL)
{
    if (!is.null(seed)) set.seed(seed)
    series_sim <- sigma_sim <- NULL
    if (nsim > 0 & h > 1) {
        if (sim_method == "bootstrap") {
            out <- garch_bootstrap(object, h = h, nsim = nsim, block = block, vreg = tail(vreg, h), seed = seed)
            sigma_sim <- out$sigma
            colnames(sigma_sim) <- as.character(forc_dates)
            class(sigma_sim) <- "tsmodel.distribution"
            attr(sigma_sim, "date_class") <- "Date"
            series_sim <- out$series
            colnames(series_sim) <- as.character(forc_dates)
            class(series_sim) <- "tsmodel.distribution"
            attr(series_sim, "date_class") <- "Date"
        } else {
            spec <- object$spec
            spec$parmatrix <- copy(object$parmatrix)
            zsim <- rdist(object$spec$distribution, h * nsim, 0, 1, skew = spec$parmatrix[parameter == "skew"]$value, shape = spec$parmatrix[parameter == "shape"]$value, lambda = spec$parmatrix[parameter == "lambda"]$value)
            zsim <- matrix(zsim, ncol = h, nrow = nsim)
            maxpq <- max(spec$model$order)
            z <- as.numeric(residuals(object, standardize = TRUE))
            init_v <- tail(as.numeric(object$sigma), maxpq)^2
            init_z <- tail(z, maxpq)
            # if model == cgarch must provide a matrix for var_init
            if (object$spec$model$model == "cgarch") {
                init_v <- cbind(init_v, tail(object$permanent_component, maxpq))
            }
            out <- simulate(spec, h = h, nsim = nsim, var_init = init_v, innov = zsim, innov_init = init_z, vreg = tail(vreg, h), seed = seed)
            sigma_sim <- out$sigma
            colnames(sigma_sim) <- as.character(forc_dates)
            class(sigma_sim) <- "tsmodel.distribution"
            attr(sigma_sim, "date_class") <- "Date"
            series_sim <- out$series
            colnames(series_sim) <- as.character(forc_dates)
            class(series_sim) <- "tsmodel.distribution"
            attr(series_sim, "date_class") <- "Date"
        }
    }
    return(list(series = series_sim, sigma = sigma_sim))
}



.predict_garch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = 1000, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = list(), seed = NULL, ...)
{
    init_model <- setup_prediction(object, h = h, newxreg = newxreg, newvreg = newvreg, forc_dates = forc_dates)
    maxpq <- init_model$maxpq
    model_parameters <- init_model$model_parameters
    constant <- model_parameters$omega + init_model$variance_regressors
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    init_states <- initialize_states(object, init_states)
    res <- c(init_states$residuals, rep(0, h))
    sigma_sqr <- c(init_states$variance, rep(0, h))
    res_sqr <- res^2
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_sqr[i] <- constant[i]
        if (init_model$garch_order[2] > 0) {
            for (j in 1:init_model$garch_order[2]) {
                sigma_sqr[i] <- sigma_sqr[i] + model_parameters$beta[j] * sigma_sqr[i - j]
            }
        }
        if (init_model$garch_order[1] > 0) {
            for (j in 1:init_model$garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- sigma_sqr[(i - j)]
                } else {
                    s <- res[(i - j)]^2
                }
                sigma_sqr[i] <- sigma_sqr[i] + model_parameters$alpha[j] * s
            }
        }
    }
    sigma <- sqrt(sigma_sqr)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    y <- y + model_parameters$mu
    series_sim <- NULL
    sigma_sim <- NULL
    simd <- simulated_distribution(object, sigma = sigma, h = h, nsim = nsim,
                                   model_parameters = model_parameters,
                                   vreg = init_model$variance_regressors,
                                   forc_dates = init_model$forc_dates,
                                   sim_method = sim_method, block = block, seed = seed)
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = simd$series,
              sigma_sim  = simd$sigma, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}


.predict_egarch_analytic <- function(object, h = 1, init_model, init_states, ...)
{
    # for 1-1 model
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    kappa <- object$kappa
    model_parameters <- init_model$model_parameters
    vreg = init_model$variance_regressors
    sigma_sqr <- rep(0, h)
    prod_approx <- .egarch_moment_prod_approx(alpha = model_parameters$alpha,
                                              gamma = model_parameters$gamma,
                                              beta = model_parameters$beta,
                                              kappa = kappa,
                                              skew = model_parameters$distribution[1],
                                              shape = model_parameters$distribution[2],
                                              lambda = model_parameters$distribution[3],
                                              distribution = object$spec$distribution,
                                              n = h)
    # initialize model:
    sigma_sqr[1] <- exp(model_parameters$omega + vreg[1]) * init_states$variance[1]^(model_parameters$beta) * exp(model_parameters$alpha * init_states$std_residuals[1] + model_parameters$gamma * (abs(init_states$std_residuals[1]) - kappa))
    if (h > 1) {
        for (i in 2:h) {
            sigma_sqr[i] <- sigma_sqr[1]^(model_parameters$beta^(i - 1)) * exp((1 - model_parameters$beta^(i - 1))/(1 - model_parameters$beta) * (model_parameters$omega + vreg[i])) * prod(prod_approx[1:(i - 1)])
        }
    }
    return(sigma_sqr)
}

.predict_egarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = NULL, seed = NULL, ...)
{
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    init_model <- setup_prediction(object, h = h, newxreg = newxreg, newvreg = newvreg, forc_dates = forc_dates)
    maxpq <- init_model$maxpq
    model_parameters <- init_model$model_parameters
    constant <- model_parameters$omega + init_model$variance_regressors
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    init_states <- initialize_states(object, init_states)
    res <- c(init_states$residuals, rep(0, h))
    sigma_sqr <- c(init_states$variance, rep(0, h))
    res_sqr <- res^2
    y <- rep(0, h)
    if (maxpq == 1) {
        sigma <- sqrt(.predict_egarch_analytic(object, h = h, init_model = init_model, init_states = init_states))
    } else {
        spec_copy <- object$spec
        spec_copy$parmatrix <- copy(object$parmatrix)
        sigma <- simulate(spec_copy, nsim = 10000, h = h, var_init = init_states$variance, innov_init = init_states$std_residuals, seed = seed)
        sigma <- sqrt(as.numeric(apply(sigma$sigma^2, 2, mean)))
    }
    y <- rep(model_parameters$mu, h)
    series_sim <- NULL
    sigma_sim <- NULL
    simd <- simulated_distribution(object, sigma = sigma, h = h, nsim = nsim, model_parameters = model_parameters,
                                   vreg = init_model$variance_regressors, forc_dates = init_model$forc_dates, sim_method = sim_method[1], block = block, seed = seed)
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = simd$series, sigma_sim  = simd$sigma, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}



.predict_aparch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = NULL, seed = NULL, ...)
{
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    init_model <- setup_prediction(object, h = h, newxreg = newxreg, newvreg = newvreg, forc_dates = forc_dates)
    maxpq <- init_model$maxpq
    model_parameters <- init_model$model_parameters
    constant <- model_parameters$omega + init_model$variance_regressors
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    init_states <- initialize_states(object, init_states)
    kappa <- object$kappa
    res <- c(init_states$residuals, rep(0, h))
    power_sigma <- c(init_states$variance^(model_parameters$delta/2), rep(0, h))

    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        power_sigma[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                power_sigma[i] <- power_sigma[i] + model_parameters$beta[j] * power_sigma[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- kappa[j] * power_sigma[i - j]
                } else {
                    s <- (abs(res[i - j]) - model_parameters$gamma[j] * res[i - j])^model_parameters$delta
                }
                power_sigma[i] <- power_sigma[i] + model_parameters$alpha[j] * s
            }
        }
    }
    sigma <- power_sigma^(1/model_parameters$delta)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    y <- y + model_parameters$mu
    series_sim <- NULL
    sigma_sim <- NULL
    simd <- simulated_distribution(object, sigma = sigma, h = h, nsim = nsim,
                                   model_parameters = model_parameters,
                                   vreg = init_model$variance_regressors,
                                   forc_dates = init_model$forc_dates,
                                   sim_method = sim_method, block = block, seed = seed)
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = simd$series,
              sigma_sim  = simd$sigma, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}

.predict_gjrgarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = NULL, seed = NULL, ...)
{
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    init_model <- setup_prediction(object, h = h, newxreg = newxreg, newvreg = newvreg, forc_dates = forc_dates)
    maxpq <- init_model$maxpq
    model_parameters <- init_model$model_parameters
    constant <- model_parameters$omega + init_model$variance_regressors
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    init_states <- initialize_states(object, init_states)
    kappa <- object$kappa
    res <- c(init_states$residuals, rep(0, h))
    sigma_sqr <- c(init_states$variance, rep(0, h))
    res_sqr <- res^2
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_sqr[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                sigma_sqr[i] <- sigma_sqr[i] + model_parameters$beta[j] * sigma_sqr[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s1 <- sigma_sqr[i - j]
                    s2 <- model_parameters$gamma[j] * kappa * sigma_sqr[i - j]
                } else {
                    s1 <- res_sqr[i - j]
                    s2 <- model_parameters$gamma[j] * (res_sqr[i - j] * (res[i - j] <= 0))
                }
                sigma_sqr[i] <- sigma_sqr[i] + model_parameters$alpha[j] * s1 + s2
            }
        }
    }
    sigma <- sqrt(sigma_sqr)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    y <- y + model_parameters$mu
    series_sim <- NULL
    sigma_sim <- NULL
    simd <- simulated_distribution(object, sigma = sigma, h = h, nsim = nsim,
                                   model_parameters = model_parameters,
                                   vreg = init_model$variance_regressors,
                                   forc_dates = init_model$forc_dates,
                                   sim_method = sim_method, block = block, seed = seed)
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = simd$series,
              sigma_sim  = simd$sigma, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}

.predict_fgarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = NULL, seed = NULL, ...)
{
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    init_model <- setup_prediction(object, h = h, newxreg = newxreg, newvreg = newvreg, forc_dates = forc_dates)
    maxpq <- init_model$maxpq
    model_parameters <- init_model$model_parameters
    constant <- model_parameters$omega + init_model$variance_regressors
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    init_states <- initialize_states(object, init_states)
    kappa <- object$kappa
    std_res <- c(init_states$std_residuals, rep(0, h))
    power_sigma <- c(init_states$variance^(model_parameters$delta/2), rep(0, h))
    y <- rep(0, h)
    for (i in (maxpq + 1):(h + maxpq)) {
        power_sigma[i] <- constant[i]
        if (garch_order[2] > 0) {
            for (j in 1:garch_order[2]) {
                power_sigma[i] <- power_sigma[i] + model_parameters$beta[j] * power_sigma[i - j]
            }
        }
        if (garch_order[1] > 0) {
            for (j in 1:garch_order[1]) {
                if ((i - maxpq) > j) {
                    s <- power_sigma[i - j] * kappa[j]
                } else {
                    s <-  power_sigma[i - j] * (abs(std_res[i - j] -  model_parameters$eta[j]) -  model_parameters$gamma[j] * (std_res[i - j] -  model_parameters$eta[j]))^model_parameters$delta
                }
                power_sigma[i] <- power_sigma[i] +  model_parameters$alpha[j] * s
            }
        }
    }
    sigma <- power_sigma^(1/model_parameters$delta)
    if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    y <- y + model_parameters$mu
    series_sim <- NULL
    sigma_sim <- NULL
    simd <- simulated_distribution(object, sigma = sigma, h = h, nsim = nsim,
                                   model_parameters = model_parameters,
                                   vreg = init_model$variance_regressors,
                                   forc_dates = init_model$forc_dates,
                                   sim_method = sim_method, block = block, seed = seed)
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = simd$series,
              sigma_sim  = simd$sigma, spec = object$spec, sigma = sigma, mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}

.predict_cgarch <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = NULL, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = NULL, seed = NULL, ...)
{
    garch_order <- object$spec$model$order
    maxpq <- max(garch_order)
    init_model <- setup_prediction(object, h = h, newxreg = newxreg, newvreg = newvreg, forc_dates = forc_dates)
    maxpq <- init_model$maxpq
    model_parameters <- init_model$model_parameters
    constant <- model_parameters$omega + init_model$variance_regressors
    if (object$spec$vreg$multiplicative) constant <- exp(constant)
    init_states <- initialize_states(object, init_states)
    res <- c(init_states$residuals, rep(0, h))
    sigma_sqr <- c(init_states$variance, rep(0, h))
    permanent_component <- c(init_states$permanent_component, rep(0, h))
    transitory_component <- c(init_states$transitory_component, rep(0, h))
    res_sqr <- res^2
    y <- rep(0, h)
    k <- 1
    for (i in (maxpq + 1):(h + maxpq)) {
        if (k == 1) {
            permanent_component[i] <- constant[i] + model_parameters$rho * permanent_component[i - 1] + model_parameters$phi * (res_sqr[i - 1] - sigma_sqr[i - 1])
            if (garch_order[1] > 0) {
                for (j in 1:garch_order[1]) {
                    if ((i - maxpq) > j) {
                        s <- sigma_sqr[i - j]
                    } else {
                        s <- res_sqr[i - j]
                    }
                    transitory_component[i] <- transitory_component[i] + model_parameters$alpha[j] * (s - sigma_sqr[i - j]) + model_parameters$alpha[j] * transitory_component[i - j]
                }
            }
            if (garch_order[2] > 0) {
                for (j in 1:garch_order[1]) {
                    transitory_component[i] <- transitory_component[i] + model_parameters$beta[j] * transitory_component[i - j]
                }
            }
            sigma_sqr[i] <- permanent_component[i] + transitory_component[i]
        } else {
            permanent_component[i] <- constant[i]/(1 - model_parameters$rho) + model_parameters$rho^k * (permanent_component[maxpq] - constant[i]/(1 - model_parameters$rho))
            p <- sum(model_parameters$beta) + sum(model_parameters$alpha)
            sigma_sqr[i] <- permanent_component[i]  +  p^k * (sigma_sqr[maxpq] - permanent_component[maxpq])
        }
        k <- k + 1
    }
    sigma <- sqrt(sigma_sqr)
    if (maxpq > 0) {
        sigma <- sigma[-seq_len(maxpq)]
        permanent_component <- permanent_component[-seq_len(maxpq)]
        transitory_component <- transitory_component[-seq_len(maxpq)]

    }

    y <- y + model_parameters$mu
    series_sim <- NULL
    sigma_sim <- NULL
    simd <- simulated_distribution(object, sigma = sigma, h = h, nsim = nsim,
                                   model_parameters = model_parameters,
                                   vreg = init_model$variance_regressors,
                                   forc_dates = init_model$forc_dates,
                                   sim_method = sim_method, block = block, seed = seed)
    mu <- xts(y, forc_dates)
    sigma <- xts(sigma, forc_dates)
    L <- list(original_series = object$spec$target$y, distribution = simd$series,
              sigma_sim  = simd$sigma, spec = object$spec, sigma = sigma,
              permanent_component = permanent_component,
              transitory_component = transitory_component,
              mean = mu)
    class(L) <- c("tsgarch.predict","tsmodel.predict")
    return(L)
}


garch_bootstrap <- function(object, h = 2, nsim = 1, block = 1, vreg = NULL, seed = seed)
{
    z <- as.numeric(residuals(object, standardize = TRUE))
    zsim <- sample_block(x = z, h = h, nsim = nsim, block = block)
    spec <- object$spec
    spec$parmatrix <- copy(object$parmatrix)
    maxpq <- max(spec$model$order)
    init_v <- tail(as.numeric(object$sigma), maxpq)^2
    init_z <- tail(z, maxpq)
    # if model == cgarch must provide a matrix for var_init
    if (object$spec$model$model == "cgarch") {
        init_v <- cbind(init_v, tail(object$permanent_component, maxpq))
    }
    b <- simulate(spec, h = h, nsim = nsim, var_init = init_v, innov = zsim, innov_init = init_z, vreg = vreg, seed = seed)
    return(b)
}

sample_block <- function(x, h, nsim, block) {
    if (block == 1) {
        zsim <- do.call(rbind, lapply(1:nsim, function(i){
            sample(x, h, replace = TRUE)
        }))
    } else {
        no_samples <- floor(h/block)
        ematrix <- .embed(as.numeric(x), block)
        n <- NROW(ematrix)
        zsim <- do.call(rbind, lapply(1:nsim, function(i){
            .sample_rows(ematrix, no_samples + 3, h = h)
        }))
    }
    # scale to avoid bias
    zsim <- scale(zsim)
    return(zsim)
}


.sample_rows <- function(x, n, h)
{
    idx <- sample(seq_len(NROW(x)), size = n, replace = TRUE)
    out <- as.vector(t(x[idx,,drop = FALSE]))
    out <- out[1:h]
    return(out)
}

.embed <- function(x, k, by = 1, ascending = TRUE)
{
    if (is.null(dim(x)[1]))
        n <- length(x)
    else n <- dim(x)[1]
    s <- seq(1, n - k + 1, by)
    lens <- length(s)
    cols <- if (ascending) 1:k else k:1
    return(matrix(x[s + rep(cols, rep(lens, k)) - 1], lens))
}


