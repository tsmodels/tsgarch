.filter.tsgarch.spec <- function(object, y = NULL, newxreg = NULL, newvreg = NULL, ...)
{
    if (!is.null(y)) {
        valid_data <- .check_y_filter(object, y = y, newvreg = newvreg)
    }
    parameter <- group <- value <- NULL
    newspec <- .spec2newspec(object, y = NULL, newxreg = NULL, newvreg = NULL)
    newspec$parmatrix$value <- copy(object$parmatrix$value)
    model <- newspec$model$model
    model_init <- .tmb_initialize_model(spec = newspec)
    tmb <- MakeADFun(data = model_init$data, parameters = model_init$parameters, atomic = TRUE, map = model_init$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    env <- new.env()
    env$fun <- model_init$fun
    env$grad <- model_init$grad
    env$hess <- model_init$hess
    env$tmb <- tmb
    env$llh <- 1
    env$model <- model
    env$distribution <- newspec$distribution
    env$parmatrix <- copy(newspec$parmatrix)
    pars <- tmb$par
    hessian <- tmb$he()
    scores <- jacobian(score_function, pars, env = env)
    m <- newspec$model_options[1]
    sig <- env$tmb$report(pars)$sigma
    if (m > 0) sig <- sig[-seq_len(m)]
    var_initial <- env$tmb$report(pars)$initial_variance
    arch_initial <- env$tmb$report(pars)$initial_arch
    target_omega <- env$tmb$report(pars)$target_omega
    rr <- suppressWarnings(summary(sdreport(tmb, par.fixed = pars, getReportCovariance = T), p.value = TRUE))
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]
    parmatrix <- copy(newspec$parmatrix)
    spec <- newspec
    spec$parmatrix <- NULL
    spec$model$var_initial <- var_initial
    constant_variance <- mean((newspec$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)
    parmatrix[parameter == "omega", value := target_omega]
    # return the ll_vector
    llvector <- -1.0 * log(tmb$report(pars)$ll_vector)
    out <- list(parmatrix = parmatrix, scaled_hessian = hessian,
                scaled_scores = scores,
                parameter_scale = rep(1, length(pars)),
                conditions = NULL,
                var_initial = var_initial,
                arch_initial = arch_initial,
                constant_variance = constant_variance,
                target_omega = target_omega,
                sigma = sig,
                loglik = tmb$fn(pars),
                lik_vector = llvector,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    if (object$model$model == "cgarch") {
        permanent_component <- env$tmb$report(pars)$permanent_component
        transitory_component <- env$tmb$report(pars)$transitory_component
        if (m > 0) {
            permanent_component <- permanent_component[-seq_len(m)]
            transitory_component <- transitory_component[-seq_len(m)]
        }
        out$permanent_component <- permanent_component
        out$transitory_component <- transitory_component
    }
    out$kappa <- NULL
    if (object$model$model %in% c("egarch","aparch","fgarch","gjrgarch")) {
        out$kappa <- env$tmb$report(pars)$kappa
    }

    class(out) <- "tsgarch.estimate"
    if (!is.null(y)) out <- tsfilter(out, y = y, newxreg = newxreg, newvreg = newvreg)
    return(out)
}

.filter_model_values <- function(object) {
    model <- object$spec$model$model
    v_orig <- extract_model_values(object, object_type = "estimate", value_name = "vreg")
    mu <- extract_model_values(object, object_type = "estimate", value_name = "mu")
    alpha <- extract_model_values(object, object_type = "estimate", value_name = "alpha")
    beta <- extract_model_values(object, object_type = "estimate", value_name = "beta")
    xi <- extract_model_values(object, object_type = "estimate", value_name = "xi")
    dpars <- extract_model_values(object, object_type = "estimate", value_name = "distribution")
    omega <- omega(object)
    L <- list(v_orig = v_orig, mu = mu, alpha = alpha, beta = beta, xi = xi, dpars = dpars, omega = omega)
    if (model == "egarch" | model == "gjrgarch") {
        gamma <- extract_model_values(object, object_type = "estimate", value_name = "gamma")
        L$gamma <- gamma
    } else if (model == "aparch") {
        gamma <- extract_model_values(object, object_type = "estimate", value_name = "gamma")
        delta <- extract_model_values(object, object_type = "estimate", value_name = "delta")
        L$gamma <- gamma
        L$delta <- delta
    } else if (model == "fgarch") {
        gamma <- extract_model_values(object, object_type = "estimate", value_name = "gamma")
        delta <- extract_model_values(object, object_type = "estimate", value_name = "delta")
        eta <- extract_model_values(object, object_type = "estimate", value_name = "eta")
        L$gamma <- gamma
        L$delta <- delta
        L$eta <- eta
    } else if (model == "cgarch") {
        rho <- extract_model_values(object, object_type = "estimate", value_name = "rho")
        phi <- extract_model_values(object, object_type = "estimate", value_name = "phi")
        L$rho <- rho
        L$phi <- phi
    } else {
        return(L)
    }
    return(L)
}

.check_y_filter <- function(object, y = NULL, newvreg = NULL)
{
    if (!is.null(y)) {
        index_new_y <- index(y)
        index_old_y <- object$target$index
        check <- max(index_old_y) < min(index_new_y)
        if (!check) {
            stop("\none of more timestamps in y is before the timestamps in the object data.")
        }
        if (!is.null(newvreg)) {
            if (object$vreg$include_vreg) {
                newvreg <- as.matrix(newvreg)
                if (NROW(newvreg) != NROW(y)) stop('\nnewvreg must have the same number of rows as y.')
            } else {
                newvreg <- NULL
            }
        } else {
            if (object$vreg$include_vreg) {
                newvreg <- as.matrix(0, ncol = ncol(object$vreg$vreg), nrow = NROW(y))
                warning('\nnewvreg is NULL but model object uses variance regressors...setting to zero.')
            } else {
                newvreg <- NULL
            }
        }
    } else {
        newvreg <- NULL
    }
    return(list(y = y, newvreg = newvreg))
}

.filter.tsgarch.estimate <- function(object, y = NULL, newxreg = NULL, newvreg = NULL, ...)
{
    # omega [init_variance, omega]
    # v is the external regressor V x \xi (already pre-multiplied in the R code)
    # model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    object_type <- "estimate"
    parameter <- group <- NULL
    if (is.null(y)) {
        return(object)
    }
    if (!is.xts(y)) stop("\ny must be an xts vector")
    if (!is.null(y)) {
        # this provides stricter checks than .merge_data
        valid_data <- .check_y_filter(object$spec, y = y, newvreg = newvreg)
        y <- valid_data$y
        newvreg <- valid_data$newvreg
    }
    maxpq <- max(object$spec$model$order)
    spec <- object$spec
    init_var <- tail(object$sigma^2, maxpq)
    y_new <- .merge_data(spec$target$y, y)
    new_n <- NROW(y_new) - NROW(spec$target$y)
    n <- NROW(y)
    maxpq <- max(spec$model$order)
    # [maxpq arch_order garch_order multiplicative]
    model <- c(maxpq, spec$model$order, as.integer(spec$vreg$multiplicative))
    L <- .filter_model_values(object)
    v_new <- .process_filter_regressors(old_regressors = L$v_orig, new_regressors = newvreg, new_index = index(y), new_n = new_n,
                                        include_regressors = spec$vreg$include_vreg)
    initstate <- init_var
    # special initialization for 2 component
    if (object$spec$model$model == "cgarch") {
        initstate <- cbind(tail(object$transitory_component, maxpq), tail(object$permanent_component, maxpq))
    }
    v <- as.numeric(v_new %*% L$xi)
    residuals <- as.numeric(y_new) - L$mu
    residuals <- tail(residuals, n + maxpq)
    v <- tail(v, n + maxpq)
    # Rcpp code
    negative_indicator <- 1 * (residuals <= 0)

    filtered_sigma <- switch(object$spec$model$model,
                    "garch"  = .garchfilter(residuals = residuals, v = v, initstate = initstate, omega = L$omega, alpha = L$alpha, beta = L$beta, model = model),
                    "egarch" = .egarchfilter(residuals = residuals, v = v, initstate = initstate, omega = L$omega, alpha = L$alpha, gamma = L$gamma, beta = L$beta,
                                             kappa = object$kappa, model = model),
                    "aparch" = .aparchfilter(residuals = residuals, v = v, initstate = initstate, omega = L$omega, alpha = L$alpha, gamma = L$gamma, beta = L$beta,
                                             delta = L$delta, model = model),
                    "fgarch" = .fgarchfilter(residuals = residuals, v = v, initstate = initstate, omega = L$omega, alpha = L$alpha, gamma = L$gamma, eta = L$eta,
                                             beta = L$beta, delta = L$delta, model = model),
                    "gjrgarch" = .gjrgarchfilter(residuals = residuals, negative_indicator = negative_indicator, v = v, initstate = initstate, omega = L$omega,
                                                 alpha = L$alpha, gamma = L$gamma, beta = L$beta, model = model),
                    "cgarch" = .cgarchfilter(residuals = residuals, v = v, initstate = initstate, omega = L$omega, alpha = L$alpha, rho = L$rho, phi = L$phi,
                                             beta = L$beta, model = model))
    if (object$spec$model$model == "cgarch") {
        sigma <- filtered_sigma$sigma
        permanent_component <- filtered_sigma$permanent_component
        transitory_component <- filtered_sigma$transitory_component
        if (maxpq > 0) {
            sigma <- sigma[-seq_len(maxpq)]
            permanent_component <- permanent_component[-seq_len(maxpq)]
            transitory_component <- transitory_component[-seq_len(maxpq)]
        }
        object$permanent_component <- c(object$permanent_component, permanent_component)
        object$transitory_component <- c(object$transitory_component, transitory_component)
    } else {
        sigma <- filtered_sigma
        if (maxpq > 0) sigma <- sigma[-seq_len(maxpq)]
    }
    object$sigma <- c(object$sigma, sigma)
    # create filter object for spec input
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    logl <- -sum(ddist(object$spec$distribution, (as.numeric(y_new) - L$mu)/object$sigma, 0, 1, skew = L$dpars[1], shape = L$dpars[2], lambda = L$dpars[3], log = TRUE) - log(object$sigma))
    object$loglik <- logl
    object$spec$target$y_orig <- as.numeric(y_new)
    # add filtered dates (increment)
    if (is.null(object$spec$target$filtered_index)) {
        object$spec$target$filtered_index <- index(y)
    } else {
        object$spec$target$filtered_index <- c(object$spec$target$filtered_index, index(y))
    }
    object$spec$target$y <- y_new
    object$spec$target$index <- index(y_new)
    object$spec$target$good <- good
    object$spec$vreg$vreg <- v_new
    object$nobs <- length(y_new)
    return(object)
}
