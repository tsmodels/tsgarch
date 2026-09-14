.filter.tsgarch.spec <- function(object, y = NULL, newxreg = NULL, newvreg = NULL, ...)
{
    if (!is.null(y)) {
        valid_data <- .check_y_filter(object, y = y, newvreg = newvreg, newxreg = newxreg)
        newvreg <- valid_data$newvreg
        newxreg <- valid_data$newxreg
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
    # conditional_mean(t), only reported by the "garch" TMB template when
    # arma > c(0,0) (see garchfun.hpp); used to seed subsequent incremental
    # tsfilter() calls (see .filter.tsgarch.estimate() / arma_filter_extend()).
    conditional_mu <- NULL
    # also capture the conditional mean when the model has mean regressors
    # with no ARMA dynamics (conditional_mean = mu + x'tau is time-varying)
    if ((!is.null(newspec$model$arma) && sum(newspec$model$arma) > 0) || isTRUE(newspec$xreg$include_xreg)) {
        conditional_mu <- tail(env$tmb$report(pars)$conditional_mean, length(sig))
    }
    # the variance target is the unconditional variance of the ARMA
    # innovations eps = y - conditional_mu over the full sample (see
    # vignettes/garch_models.Rmd, "Variance Targeting"); fall back to the
    # constant-mean deviation when no ARMA mean equation is present
    if (!is.null(conditional_mu)) {
        constant_variance <- mean((as.numeric(newspec$target$y_orig) - conditional_mu)^2)
    } else {
        constant_variance <- mean((newspec$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)
    }
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
                spec = spec,
                conditional_mu = conditional_mu)
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
    tau <- extract_model_values(object, object_type = "estimate", value_name = "tau")
    x_orig <- extract_model_values(object, object_type = "estimate", value_name = "xreg")
    # older objects may lack tau rows entirely; pad to the regressor column
    # count so the matrix product below is always conformable
    if (length(tau) != NCOL(x_orig)) tau <- rep(0, NCOL(x_orig))
    armax <- isTRUE(object$spec$xreg$xreg_type == "armax")
    dpars <- extract_model_values(object, object_type = "estimate", value_name = "distribution")
    omega <- omega(object)
    arma_order <- object$spec$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    if (sum(arma_order) > 0) {
        arma_coef <- arma_coefficients(object)
        ar <- as.numeric(arma_coef$ar)
        ma <- as.numeric(arma_coef$ma)
    } else {
        ar <- numeric(0)
        ma <- numeric(0)
    }
    L <- list(v_orig = v_orig, mu = mu, alpha = alpha, beta = beta, xi = xi, tau = tau, x_orig = x_orig, armax = armax, dpars = dpars, omega = omega, ar = ar, ma = ma)
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

.check_y_filter <- function(object, y = NULL, newvreg = NULL, newxreg = NULL)
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
                if (any(!is.finite(newvreg))) stop('\nNA/NaN/Inf values found in newvreg.')
            } else {
                newvreg <- NULL
            }
        } else {
            if (object$vreg$include_vreg) {
                newvreg <- matrix(0, nrow = NROW(y), ncol = ncol(object$vreg$vreg))
                warning('\nnewvreg is NULL but model object uses variance regressors...setting to zero.')
            } else {
                newvreg <- NULL
            }
        }
        # mean equation regressors follow the same warn-and-zero policy;
        # isTRUE() keeps the check working on objects serialized by package
        # versions predating the xreg slot
        include_xreg <- isTRUE(object$xreg$include_xreg)
        if (!is.null(newxreg)) {
            if (include_xreg) {
                newxreg <- as.matrix(newxreg)
                if (NROW(newxreg) != NROW(y)) stop('\nnewxreg must have the same number of rows as y.')
                if (any(!is.finite(newxreg))) stop('\nNA/NaN/Inf values found in newxreg.')
            } else {
                newxreg <- NULL
            }
        } else {
            if (include_xreg) {
                newxreg <- matrix(0, nrow = NROW(y), ncol = ncol(object$xreg$xreg))
                warning('\nnewxreg is NULL but model object uses mean regressors...setting to zero.')
            } else {
                newxreg <- NULL
            }
        }
    } else {
        newvreg <- NULL
        newxreg <- NULL
    }
    return(list(y = y, newvreg = newvreg, newxreg = newxreg))
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
        valid_data <- .check_y_filter(object$spec, y = y, newvreg = newvreg, newxreg = newxreg)
        y <- valid_data$y
        newvreg <- valid_data$newvreg
        newxreg <- valid_data$newxreg
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
                                        include_regressors = spec$vreg$include_vreg, regressor_argument = "newvreg")
    x_new <- .process_filter_regressors(old_regressors = L$x_orig, new_regressors = newxreg, new_index = index(y), new_n = new_n,
                                        include_regressors = isTRUE(spec$xreg$include_xreg), regressor_argument = "newxreg")
    # per-period mean regressor contribution over the full merged series
    xtau_full <- as.numeric(x_new %*% L$tau)
    initstate <- init_var
    # special initialization for 2 component
    if (object$spec$model$model == "cgarch") {
        initstate <- cbind(tail(object$transitory_component, maxpq), tail(object$permanent_component, maxpq))
    }
    v <- as.numeric(v_new %*% L$xi)
    arma_order <- spec$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    if (sum(arma_order) > 0) {
        # continue the ARMA mean recursion (see garchfun.hpp) for the newly
        # appended observations only, using the already-computed historical
        # conditional mean as the continuation state (object$conditional_mu,
        # set by estimate()/.filter.tsgarch.spec() and updated below so
        # repeated tsfilter() calls keep chaining correctly).
        old_conditional_mu <- object$conditional_mu
        if (is.null(old_conditional_mu)) {
            stop("\nobject does not carry conditional_mu but has arma > c(0,0); re-fit/re-filter the base object with the updated tsgarch version.")
        }
        new_conditional_mu <- arma_filter_extend(as.numeric(y_new), length(old_conditional_mu), L$mu, L$ar, L$ma, old_conditional_mu,
                                                 xtau_full = xtau_full, armax = L$armax)
        full_conditional_mu <- c(old_conditional_mu, new_conditional_mu)
    } else {
        # regression with GARCH errors and no ARMA dynamics still has a
        # time-varying conditional mean through x'tau
        full_conditional_mu <- rep(L$mu, NROW(y_new)) + xtau_full
    }
    full_residuals <- as.numeric(y_new) - full_conditional_mu
    residuals <- tail(full_residuals, n + maxpq)
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
    if (sum(arma_order) > 0 || isTRUE(spec$xreg$include_xreg)) object$conditional_mu <- full_conditional_mu
    # create filter object for spec input
    good <- rep(1, NROW(y_new))
    if (any(is.na(y_new))) {
        good[which(is.na(y_new))] <- 0
    }
    logl <- -sum(ddist(object$spec$distribution, full_residuals/object$sigma, 0, 1, skew = L$dpars[1], shape = L$dpars[2], lambda = L$dpars[3], log = TRUE) - log(object$sigma))
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
    object$spec$xreg$xreg <- x_new
    object$nobs <- length(y_new)
    return(object)
}
