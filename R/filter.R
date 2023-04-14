.filter_garch <- function(object, y = NULL, newxreg = NULL, newvreg = NULL, object_type = "estimate", ...)
{
    # omega [init_variance, omega]
    # v is the external regressor V x \xi (already pre-multiplied in the R code)
    # model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    parameter <- group <- NULL
    if (is.null(y)) {
        if (object_type == "estimate") {
            return(object)
        } else {
            y <- object$target$y
            newvreg <- object$vreg$vreg
        }
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (object_type == "estimate") {
        spec <- object$spec
        init_var <- object$var_initial
    } else {
        spec <- object
        init_var <- initialize_variance(y = spec$target$y_orig, mu = spec$parmatrix[parameter == "mu"]$value,
                                        init = spec$model$init, backcast_lambda = spec$model$backcast_lambda,
                                        sample_n = spec$model$sample_n)
    }
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    model <- c(maxpq, as.integer(spec$vreg$multiplicative), spec$model$order)
    v_orig <- extract_model_values(object, object_type = object_type, value_name = "vreg")
    mu <- extract_model_values(object, object_type = object_type, value_name = "mu")
    alpha <- extract_model_values(object, object_type = object_type, value_name = "alpha")
    beta <- extract_model_values(object, object_type = object_type, value_name = "beta")
    xi <- extract_model_values(object, object_type = object_type, value_name = "xi")
    dpars <- extract_model_values(object, object_type = object_type, value_name = "distribution")

    v_new <- .process_filter_regressors(old_regressors = v_orig, new_regressors = newvreg, new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    omega <- omega(object)
    # create var_initial functions in R for spec object
    initstate <- rep(init_var, maxpq)
    v <- as.numeric(v_new %*% xi)
    residuals <- as.numeric(y_new) - mu
    residuals <- c(rep(0, maxpq), residuals)
    v <- c(rep(0, maxpq), v)
    # Rcpp code
    sigma <- .garchfilter(residuals, v, initstate, omega, alpha, beta, model)
    sigma <- sigma[-c(1:maxpq)]
    # create filter object for spec input
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    if (object_type == "estimate") {
        object$sigma <- sigma
        object$spec$target$y_orig <- as.numeric(y_new)
        object$spec$target$y <- y_new
        object$spec$target$index <- index(y_new)
        object$spec$target$good <- good
        object$spec$vreg$vreg <- v_new
        object$nobs <- length(residuals) - maxpq
        return(object)
    } else {
        spec$target$y_orig <- as.numeric(y_new)
        spec$target$index <- index(y_new)
        spec$target$y <- y_new
        spec$target$good <- good
        spec$vreg$vreg <- v_new
        L <- list()
        L$parmatrix <- spec$parmatrix
        L$scaled_hessian <- NULL
        L$scaled_scores <- NULL
        L$parameter_scale <- L$parmatrix[estimate == 1]$scale
        L$conditions <- NULL
        L$var_initial <- init_var
        L$constant_variance <- mean((L$target$y_orig - (L$parmatrix[parameter == "mu"]$value * L$parmatrix[parameter == "mu"]$scale))^2)
        L$target_omega <- omega
        L$sigma <- sigma
        L$loglik <- -1 * likelihood_fun(spec$target$y_orig, distribution = spec$distribution, mu = mu, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3])
        L$nobs <- length(y_new)
        L$persistence_summary <- NULL
        L$variance_target_summary <- NULL
        L$npars <- NROW(L$parmatrix[estimate == 1]) + 1
        L$spec <- spec
        L$spec$parmatrix <- NULL
        class(L) <- "tsgarch.estimate"
        return(L)
    }
}


.filter_egarch <- function(object, y, newxreg = NULL, newvreg = NULL, object_type = "estimate", ...)
{
    parameter <- group <- NULL
    if (is.null(y)) {
        if (object_type == "estimate") {
            return(object)
        } else {
            y <- object$target$y
            newvreg <- object$vreg$vreg
        }
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (object_type == "estimate") {
        spec <- object$spec
        init_var <- object$var_initial
        kappa <- object$kappa
    } else {
        spec <- object
        init_var <- initialize_variance(y = spec$target$y_orig, mu = spec$parmatrix[parameter == "mu"]$value,
                                        init = spec$model$init, backcast_lambda = spec$model$backcast_lambda,
                                        sample_n = spec$model$sample_n)
        kappa <- egarch_moment(spec$distribution, skew = spec$parmatrix[parameter == "skew"]$value,
                               shape = spec$parmatrix[parameter == "shape"]$value,
                               lambda = spec$parmatrix[parameter == "lambda"]$value)
    }
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    model <- c(maxpq, as.integer(spec$vreg$multiplicative), spec$model$order)
    v_orig <- extract_model_values(spec, object_type = object_type, value_name = "vreg")
    mu <- extract_model_values(spec, object_type = object_type, value_name = "mu")
    alpha <- extract_model_values(spec, object_type = object_type, value_name = "alpha")
    gamma <- extract_model_values(spec, object_type = object_type, value_name = "gamma")
    beta <- extract_model_values(spec, object_type = object_type, value_name = "beta")
    xi <- extract_model_values(spec, object_type = object_type, value_name = "xi")
    dpars <- extract_model_values(object, object_type = object_type, value_name = "distribution")
    v_new <- .process_filter_regressors(old_regressors = v_orig, new_regressors = newvreg, new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    omega <- omega(object)
    initstate <- rep(init_var, maxpq)
    v <- as.numeric(v_new %*% xi)
    residuals <- as.numeric(y_new) - mu
    residuals <- c(rep(0, maxpq), residuals)
    v <- c(rep(0, maxpq), v)
    # Rcpp code
    sigma <- .egarchfilter(residuals = residuals, v = v, initstate = initstate, omega = omega,
                           alpha = alpha, gamma = gamma, beta = beta,
                           kappa = kappa, model = model)

    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    if (object_type == "estimate") {
        object$sigma <- sigma
        object$spec$target$y_orig <- as.numeric(y_new)
        object$spec$target$y <- y_new
        object$spec$target$index <- index(y_new)
        object$spec$target$good <- good
        object$spec$vreg$vreg <- v_new
        object$nobs <- length(residuals) - maxpq
        return(object)
    } else {
        spec$target$y_orig <- as.numeric(y_new)
        spec$target$index <- index(y_new)
        spec$target$y <- y_new
        spec$target$good <- good
        spec$vreg$vreg <- v_new
        L <- list()
        L$parmatrix <- spec$parmatrix
        L$scaled_hessian <- NULL
        L$scaled_scores <- NULL
        L$parameter_scale <- L$parmatrix[estimate == 1]$scale
        L$conditions <- NULL
        L$var_initial <- init_var
        L$constant_variance <- mean((L$target$y_orig - (L$parmatrix[parameter == "mu"]$value * L$parmatrix[parameter == "mu"]$scale))^2)
        L$target_omega <- omega
        L$kappa <- kappa
        L$sigma <- sigma
        L$loglik <- -1 * likelihood_fun(spec$target$y_orig, distribution = spec$distribution, mu = mu, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3])
        L$nobs <- length(y_new)
        L$persistence_summary <- NULL
        L$variance_target_summary <- NULL
        L$kappa_summary <- NULL
        L$npars <- NROW(L$parmatrix[estimate == 1]) + 1
        L$spec <- spec
        L$spec$parmatrix <- NULL
        class(L) <- "tsgarch.estimate"
        return(L)
    }
}


.filter_aparch <- function(object, y, newxreg = NULL, newvreg = NULL, object_type = "estimate", ...)
{
    parameter <- group <- NULL
    if (is.null(y)) {
        if (object_type == "estimate") {
            return(object)
        } else {
            y <- object$target$y
            newvreg <- object$vreg$vreg
        }
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (object_type == "estimate") {
        spec <- object$spec
        init_var <- object$var_initial
        kappa <- object$kappa
    } else {
        spec <- object
        init_var <- initialize_variance(y = spec$target$y_orig, mu = spec$parmatrix[parameter == "mu"]$value,
                                        init = spec$model$init, backcast_lambda = spec$model$backcast_lambda,
                                        sample_n = spec$model$sample_n, delta = spec$parmatrix[parameter == "delta"]$value)
        kappa <- aparch_moment_v(distribution = spec$distribution, gamma = spec$parmatrix[group == "gamma"]$value,
                                 delta = spec$parmatrix[group == "delta"]$value,
                                 skew = spec$parmatrix[parameter == "skew"]$value,
                                 shape = spec$parmatrix[parameter == "shape"]$value,
                                 lambda = spec$parmatrix[parameter == "lambda"]$value)
    }
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    model <- c(maxpq, as.integer(spec$vreg$multiplicative), spec$model$order)

    v_orig <- extract_model_values(object, object_type = object_type, value_name = "vreg")
    mu <- extract_model_values(object, object_type = object_type, value_name = "mu")
    alpha <- extract_model_values(object, object_type = object_type, value_name = "alpha")
    gamma <- extract_model_values(object, object_type = object_type, value_name = "gamma")
    beta <- extract_model_values(object, object_type = object_type, value_name = "beta")
    delta <- extract_model_values(object, object_type = object_type, value_name = "delta")
    xi <- extract_model_values(object, object_type = object_type, value_name = "xi")
    dpars <- extract_model_values(object, object_type = object_type, value_name = "distribution")
    v_new <- .process_filter_regressors(old_regressors = v_orig, new_regressors = newvreg,
                                        new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    omega <- omega(object)
    initstate <- rep(init_var, maxpq)
    v <- as.numeric(v_new %*% xi)
    residuals <- as.numeric(y_new) - mu
    residuals <- c(rep(0, maxpq), residuals)
    v <- c(rep(0, maxpq), v)
    # Rcpp code
    sigma <- .aparchfilter(residuals = residuals, v = v, initstate = initstate,
                           omega = omega, alpha = alpha, gamma = gamma, beta = beta,
                           delta = delta, model = model)
    sigma <- sigma[-c(1:maxpq)]
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    if (object_type == "estimate") {
        object$sigma <- sigma
        object$spec$target$y_orig <- as.numeric(y_new)
        object$spec$target$y <- y_new
        object$spec$target$index <- index(y_new)
        object$spec$target$good <- good
        object$spec$vreg$vreg <- v_new
        object$nobs <- length(residuals) - maxpq
        return(object)
    } else {
        spec$target$y_orig <- as.numeric(y_new)
        spec$target$index <- index(y_new)
        spec$target$y <- y_new
        spec$target$good <- good
        spec$vreg$vreg <- v_new
        L <- list()
        L$parmatrix <- spec$parmatrix
        L$scaled_hessian <- NULL
        L$scaled_scores <- NULL
        L$parameter_scale <- L$parmatrix[estimate == 1]$scale
        L$conditions <- NULL
        L$var_initial <- init_var
        L$constant_variance <- mean((L$target$y_orig - (L$parmatrix[parameter == "mu"]$value * L$parmatrix[parameter == "mu"]$scale))^2)
        L$target_omega <- omega
        L$kappa <- kappa
        L$sigma <- sigma
        L$loglik <- -1 * likelihood_fun(spec$target$y_orig, distribution = spec$distribution, mu = mu, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3])
        L$nobs <- length(y_new)
        L$persistence_summary <- NULL
        L$variance_target_summary <- NULL
        L$kappa_summary <- NULL
        L$npars <- NROW(L$parmatrix[estimate == 1]) + 1
        L$spec <- spec
        L$spec$parmatrix <- NULL
        class(L) <- "tsgarch.estimate"
        return(L)
    }
}

.filter_gjrgarch <- function(object, y, newxreg = NULL, newvreg = NULL, object_type = "estimate", ...)
{
    parameter <- group <- NULL
    if (is.null(y)) {
        if (object_type == "estimate") {
            return(object)
        } else {
            y <- object$target$y
            newvreg <- object$vreg$vreg
        }
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (object_type == "estimate") {
        spec <- object$spec
        init_var <- object$var_initial
        kappa <- object$kappa
    } else {
        spec <- object
        init_var <- initialize_variance(y = spec$target$y_orig, mu = spec$parmatrix[parameter == "mu"]$value,
                                        init = spec$model$init, backcast_lambda = spec$model$backcast_lambda,
                                        sample_n = spec$model$sample_n, delta = 2)
        kappa <- gjrgarch_moment(distribution = spec$distribution, skew = spec$parmatrix[parameter == "skew"]$value,
                                 shape = spec$parmatrix[parameter == "shape"]$value,
                                 lambda = spec$parmatrix[parameter == "lambda"]$value)
    }
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    model <- c(maxpq, as.integer(spec$vreg$multiplicative), spec$model$order)

    v_orig <- extract_model_values(spec, object_type = object_type, value_name = "vreg")
    mu <- extract_model_values(spec, object_type = object_type, value_name = "mu")
    alpha <- extract_model_values(spec, object_type = object_type, value_name = "alpha")
    gamma <- extract_model_values(spec, object_type = object_type, value_name = "gamma")
    beta <- extract_model_values(spec, object_type = object_type, value_name = "beta")
    xi <- extract_model_values(spec, object_type = object_type, value_name = "xi")
    dpars <- extract_model_values(object, object_type = object_type, value_name = "distribution")

    v_new <- .process_filter_regressors(old_regressors = v_orig, new_regressors = newvreg, new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    omega <- omega(object)
    initstate <- rep(init_var, maxpq)
    v <- as.numeric(v_new %*% xi)
    residuals <- as.numeric(y_new) - mu
    residuals <- c(rep(0, maxpq), residuals)
    negative_indicator <- 1 * (residuals <= 0)
    if (maxpq > 0) negative_indicator[1:maxpq] <- 1
    v <- c(rep(0, maxpq), v)
    # Rcpp code
    sigma <- .gjrgarchfilter(residuals = residuals, negative_indicator, v = v, initstate = initstate,
                             omega = omega, alpha = alpha, gamma = gamma, beta = beta, model = model)
    sigma <- sigma[-c(1:maxpq)]
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    if (object_type == "estimate") {
        object$sigma <- sigma
        object$spec$target$y_orig <- as.numeric(y_new)
        object$spec$target$y <- y_new
        object$spec$target$index <- index(y_new)
        object$spec$target$good <- good
        object$spec$vreg$vreg <- v_new
        object$nobs <- length(residuals) - maxpq
        return(object)
    } else {
        spec$target$y_orig <- as.numeric(y_new)
        spec$target$index <- index(y_new)
        spec$target$y <- y_new
        spec$target$good <- good
        spec$vreg$vreg <- v_new
        L <- list()
        L$parmatrix <- spec$parmatrix
        L$scaled_hessian <- NULL
        L$scaled_scores <- NULL
        L$parameter_scale <- L$parmatrix[estimate == 1]$scale
        L$conditions <- NULL
        L$var_initial <- init_var
        L$constant_variance <- mean((L$target$y_orig - (L$parmatrix[parameter == "mu"]$value * L$parmatrix[parameter == "mu"]$scale))^2)
        L$target_omega <- omega
        L$kappa <- kappa
        L$sigma <- sigma
        L$loglik <- -1 * likelihood_fun(spec$target$y_orig, distribution = spec$distribution, mu = mu, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3])
        L$nobs <- length(y_new)
        L$persistence_summary <- NULL
        L$variance_target_summary <- NULL
        L$kappa_summary <- NULL
        L$npars <- NROW(L$parmatrix[estimate == 1]) + 1
        L$spec <- spec
        L$spec$parmatrix <- NULL
        class(L) <- "tsgarch.estimate"
        return(L)
    }
    return(object)
}


.filter_fgarch <- function(object, y, newxreg = NULL, newvreg = NULL, object_type = "estimate", ...)
{
    parameter <- group <- NULL
    if (is.null(y)) {
        if (object_type == "estimate") {
            return(object)
        } else {
            y <- object$target$y
            newvreg <- object$vreg$vreg
        }
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (object_type == "estimate") {
        spec <- object$spec
        init_var <- object$var_initial
        kappa <- object$kappa
    } else {
        spec <- object
        init_var <- initialize_variance(y = spec$target$y_orig, mu = spec$parmatrix[parameter == "mu"]$value,
                                        init = spec$model$init, backcast_lambda = spec$model$backcast_lambda,
                                        sample_n = spec$model$sample_n, delta = spec$parmatrix[parameter == "delta"]$value)
        kappa <- fgarch_moment_v(distribution = spec$distribution, gamma = spec$parmatrix[group == "gamma"]$value,
                                 delta = spec$parmatrix[group == "delta"]$value,
                                 eta = spec$parmatrix[group == "eta"]$value,
                                 skew = spec$parmatrix[parameter == "skew"]$value,
                                 shape = spec$parmatrix[parameter == "shape"]$value,
                                 lambda = spec$parmatrix[parameter == "lambda"]$value)
    }
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    model <- c(maxpq, as.integer(spec$vreg$multiplicative), spec$model$order)

    v_orig <- extract_model_values(spec, object_type = object_type, value_name = "vreg")
    mu <- extract_model_values(spec, object_type = object_type, value_name = "mu")
    alpha <- extract_model_values(spec, object_type = object_type, value_name = "alpha")
    gamma <- extract_model_values(spec, object_type = object_type, value_name = "gamma")
    eta <- extract_model_values(spec, object_type = object_type, value_name = "eta")
    beta <- extract_model_values(spec, object_type = object_type, value_name = "beta")
    delta <- extract_model_values(spec, object_type = object_type, value_name = "delta")
    xi <- extract_model_values(spec, object_type = object_type, value_name = "xi")
    dpars <- extract_model_values(object, object_type = object_type, value_name = "distribution")

    v_new <- .process_filter_regressors(old_regressors = v_orig, new_regressors = newvreg, new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    omega <- omega(object)
    initstate <- rep(init_var, maxpq)
    v <- as.numeric(v_new %*% xi)
    residuals <- as.numeric(y_new) - mu
    residuals <- c(rep(0, maxpq), residuals)
    v <- c(rep(0, maxpq), v)
    # Rcpp code
    sigma <- .fgarchfilter(residuals = residuals, v = v, initstate = initstate, omega = omega,
                           alpha = alpha, gamma = gamma, eta = eta, beta = beta,
                           delta = delta, model = model)
    sigma <- sigma[-c(1:maxpq)]
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    if (object_type == "estimate") {
        object$sigma <- sigma
        object$spec$target$y_orig <- as.numeric(y_new)
        object$spec$target$y <- y_new
        object$spec$target$index <- index(y_new)
        object$spec$target$good <- good
        object$spec$vreg$vreg <- v_new
        object$nobs <- length(residuals) - maxpq
        return(object)
    } else {
        spec$target$y_orig <- as.numeric(y_new)
        spec$target$index <- index(y_new)
        spec$target$y <- y_new
        spec$target$good <- good
        spec$vreg$vreg <- v_new
        L <- list()
        L$parmatrix <- spec$parmatrix
        L$scaled_hessian <- NULL
        L$scaled_scores <- NULL
        L$parameter_scale <- L$parmatrix[estimate == 1]$scale
        L$conditions <- NULL
        L$var_initial <- init_var
        L$constant_variance <- mean((L$target$y_orig - (L$parmatrix[parameter == "mu"]$value * L$parmatrix[parameter == "mu"]$scale))^2)
        L$target_omega <- omega
        L$kappa <- kappa
        L$sigma <- sigma
        L$loglik <- -1 * likelihood_fun(spec$target$y_orig, distribution = spec$distribution, mu = mu, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3])
        L$nobs <- length(y_new)
        L$persistence_summary <- NULL
        L$variance_target_summary <- NULL
        L$kappa_summary <- NULL
        L$npars <- NROW(L$parmatrix[estimate == 1]) + 1
        L$spec <- spec
        L$spec$parmatrix <- NULL
        class(L) <- "tsgarch.estimate"
        return(L)
    }
}

.filter_cgarch <- function(object, y, newxreg = NULL, newvreg = NULL, object_type = "estimate", ...)
{
    # omega [init_variance, omega]
    # v is the external regressor V x \xi (already pre-multiplied in the R code)
    # model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    parameter <- group <- NULL
    if (is.null(y)) {
        if (object_type == "estimate") {
            return(object)
        } else {
            y <- object$target$y
            newvreg <- object$vreg$vreg
        }
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (object_type == "estimate") {
        spec <- object$spec
        init_var <- object$var_initial
    } else {
        spec <- object
        init_var <- initialize_variance(y = spec$target$y_orig, mu = spec$parmatrix[parameter == "mu"]$value,
                                        init = spec$model$init, backcast_lambda = spec$model$backcast_lambda,
                                        sample_n = spec$model$sample_n, delta = 2)
    }
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    model <- c(maxpq, as.integer(spec$vreg$multiplicative), spec$model$order)

    v_orig <- extract_model_values(spec, object_type = object_type, value_name = "vreg")
    mu <- extract_model_values(spec, object_type = object_type, value_name = "mu")
    rho <- extract_model_values(spec, object_type = object_type, value_name = "rho")
    phi <- extract_model_values(spec, object_type = object_type, value_name = "phi")
    alpha <- extract_model_values(spec, object_type = object_type, value_name = "alpha")
    beta <- extract_model_values(spec, object_type = object_type, value_name = "beta")
    xi <- extract_model_values(spec, object_type = object_type, value_name = "xi")
    dpars <- extract_model_values(object, object_type = object_type, value_name = "distribution")

    v_new <- .process_filter_regressors(old_regressors = v_orig, new_regressors = newvreg, new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    omega <- omega(object)
    initstate <- as.matrix(cbind(rep(init_var, maxpq), rep(omega/(1.0 - rho), maxpq)))
    v <- as.numeric(v_new %*% xi)
    residuals <- as.numeric(y_new) - mu
    residuals <- c(rep(0, maxpq), residuals)
    v <- c(rep(0, maxpq), v)
    # Rcpp code
    out <- .cgarchfilter(residuals = residuals, v = v, initstate = initstate,
                          omega = omega, alpha = alpha, rho = rho, phi = phi,
                          beta = beta, model = as.integer(model))
    sigma <- out$sigma
    permanent_component <- out$permanent_component
    sigma <- sigma[-c(1:maxpq)]
    permanent_component <- permanent_component[-c(1:maxpq)]
    transitory_component <- sigma - permanent_component
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    if (object_type == "estimate") {
        object$sigma <- sigma
        object$permanent_component <- permanent_component
        object$transitory_component <- transitory_component
        object$spec$target$y_orig <- as.numeric(y_new)
        object$spec$target$y <- y_new
        object$spec$target$index <- index(y_new)
        object$spec$vreg$vreg <- v_new
        object$nobs <- length(residuals) - maxpq
        return(object)
    } else {
        spec$target$y_orig <- as.numeric(y_new)
        spec$target$index <- index(y_new)
        spec$target$y <- y_new
        spec$target$good <- good
        spec$vreg$vreg <- v_new
        L <- list()
        L$parmatrix <- spec$parmatrix
        L$scaled_hessian <- NULL
        L$scaled_scores <- NULL
        L$parameter_scale <- L$parmatrix[estimate == 1]$scale
        L$conditions <- NULL
        L$var_initial <- init_var
        L$constant_variance <- mean((L$target$y_orig - (L$parmatrix[parameter == "mu"]$value * L$parmatrix[parameter == "mu"]$scale))^2)
        L$target_omega <- omega
        L$kappa <- kappa
        L$sigma <- sigma
        L$permanent_component <- permanent_component
        L$transitory_component <- transitory_component
        L$loglik <- -1 * likelihood_fun(spec$target$y_orig, distribution = spec$distribution, mu = mu, sigma = sigma, skew = dpars[1], shape = dpars[2], lambda = dpars[3])
        L$nobs <- length(y_new)
        L$persistence_summary <- NULL
        L$variance_target_summary <- NULL
        L$kappa_summary <- NULL
        L$npars <- NROW(L$parmatrix[estimate == 1]) + 1
        L$spec <- spec
        L$spec$parmatrix <- NULL
        class(L) <- "tsgarch.estimate"
        return(L)
    }
}


