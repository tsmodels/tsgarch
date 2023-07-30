#--------------------------------
# common wrappers for the objective, gradient and hessian functions
garch_fun <- function(pars, env)
{
    parameter <- value <- NULL
    parmatrix <- env$parmatrix
    parmatrix[estimate == 1, value := pars]
    llh <- env$tmb$fn(pars)
    if (!is.finite(llh) | is.na(llh)) {
        llh <- env$llh + 0.1 * abs(env$llh)
        env$llh <- llh
        return(llh)
    } else {
        env$llh <- llh
        return(llh)
    }
}

garch_grad <- function(pars, env)
{
    g <- env$tmb$gr(pars)
    if (any(is.na(g))) {
        g[1,which(is.na(g))] <- 1e10
    }
    return(g)
}

garch_hess <- function(pars, env)
{
    env$tmb$he(pars)
}

#--------------------------------
# garch model
.tmb_initializa_garch <- function(spec, ...)
{
    include <- value <- group <- .N <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    active_constraints_index <- parmatrix[estimate == 1]
    active_constraints_index <- which(active_constraints_index$group %in% c("alpha","beta"))
    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda, samplen = spec$model$sample_n, initmethod = spec$model$init, pscale = pscale, cmodel = cmodel, model = "garch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map, active_constraints_index = active_constraints_index))
}


.estimate_garch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- group <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_garch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, atomic = TRUE, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "garch"
    env$stationarity_constraint <- stationarity_constraint
    env$active_constraints_index <- L$active_constraints_index
    env$parmatrix <- object$parmatrix
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale
    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad, eval_g_ineq = garch_ineq, eval_jac_g_ineq = garch_ineq_jac, lb = lower, ub = upper, env = env, opts = control)
    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    flag <- FALSE
    if (any(pmatrix[estimate == 1 & group == "alpha"]$value == pmatrix[estimate == 1 & group == "alpha"]$lower)) {
        flag <- TRUE
        pmatrix[estimate == 1 & group == "alpha" & value == lower, value := 0.01]
        pmatrix[estimate == 1 & group == "alpha" & value == lower, lower := lower + 1e-6]
    }
    if (any(pmatrix[estimate == 1 & group == "beta"]$value == pmatrix[estimate == 1 & group == "beta"]$lower)) {
        flag <- TRUE
        pmatrix[estimate == 1 & group == "beta" & value == lower, value := 0.8]
        pmatrix[estimate == 1 & group == "beta" & value == lower, lower := lower + 1e-6]
    }
    # special case trap
    if (flag) {
        warning("\nfound potential problem in optimization with at least one parameter stuck in lower bound. Adjusting bounds and re-estimating to check.")
        lower <- pmatrix[estimate == 1]$lower * pmatrix[estimate == 1]$scale
        upper <- pmatrix[estimate == 1]$upper * 1/pmatrix[estimate == 1]$scale
        sol_check <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad, eval_g_ineq = garch_ineq, eval_jac_g_ineq = garch_ineq_jac, lb = lower, ub = upper, env = env, opts = control)
        if (sol_check$objective < sol$objective) {
            sol <- sol_check
            env$parmatrix <- pmatrix
            object$parmatrix <- pmatrix
        }
    }
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_garch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)
    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table
    parmatrix[parameter == "omega", value := scaled_solution$target_omega]
    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                sigma = scaled_solution$sigma,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    return(out)
}

.estimate_garch_scaled <- function(pars, H, object, control, stationarity_constraint)
{
    parameter <- value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]

    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_L <- .tmb_initializa_garch(scaled_object)
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "garch"
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$active_constraints_index <- scaled_L$active_constraints_index
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = garch_ineq, eval_jac_g_ineq = garch_ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he()
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) sig <- sig[-seq_len(m)]
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega
    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]

    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega,
                persistence_table = persistence_table,
                variance_target_table = variance_target_table,
                var_initial = var_initial, sigma = sig))
}

#--------------------------------
# egarch model
.tmb_initializa_egarch <- function(spec, ...)
{
    include <- parameter <- group <- value <- .N <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    active_constraints_index <- parmatrix[estimate == 1]
    active_constraints_index <- which(active_constraints_index$group %in% c("beta"))
    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda, samplen = spec$model$sample_n, initmethod = spec$model$init, pscale = pscale, cmodel = cmodel, model = "egarch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map, active_constraints_index = active_constraints_index))
}

.estimate_egarch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_egarch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "egarch"
    env$stationarity_constraint <- stationarity_constraint
    env$active_constraints_index <- L$active_constraints_index
    env$parmatrix <- object$parmatrix
    env$distribution <- object$distribution
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale

    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad,
                  eval_g_ineq = garch_ineq, eval_jac_g_ineq = egarch_ineq_jac,
                  lb = lower, ub = upper, env = env, opts = control)

    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    # we don't need special flags for egarch
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_egarch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))

    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)

    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table
    kappa_table <- scaled_solution$kappa_table
    parmatrix[parameter == "omega", value := scaled_solution$target_omega]

    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                kappa = scaled_solution$kappa,
                sigma = scaled_solution$sigma,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                kappa_summary = kappa_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    rm(tmb)
    return(out)
}

.estimate_egarch_scaled <- function(pars, H, object, control, stationarity_constraint = 0.999)
{
    value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H  <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]
    scaled_L <- .tmb_initializa_egarch(scaled_object)
    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale

    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, atomic = TRUE, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "egarch"
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_env$active_constraints_index <- scaled_L$active_constraints_index
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = garch_ineq, eval_jac_g_ineq = egarch_ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he(scaled_sol$solution)
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) sig <- sig[-seq_len(m)]
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega
    kappa <- scaled_env$tmb$report(scaled_sol$solution)$kappa

    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]
    kappa_table <- rr["kappa",]

    rm(scaled_tmb)
    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega, kappa = kappa,
                persistence_table = persistence_table,
                variance_target_table = variance_target_table,
                kappa_table = kappa_table,
                var_initial = var_initial, sigma = sig))
}

#--------------------------------
# aparch model
aparch_jac_fun <- function(L, parmatrix)
{
    C <- list()
    C$data$cmodel <- L$data$cmodel
    C$data$model <- "aparchjacobian"
    C$parameters <- L$parameters
    C$parameters$mu <- NULL
    C$parameters$omega <- NULL
    C$parameters$xi <- NULL
    C$map <- L$map
    C$map$mu <- NULL
    C$map$omega <- NULL
    C$map$xi <- NULL
    n <- sum(parmatrix$estimate)
    jacobian_matrix <- matrix(0, ncol = sum(parmatrix$estimate), nrow = 1)
    groups <- parmatrix$group
    active_constraints_index <- which(parmatrix[estimate == 1]$group %in% c("alpha","gamma","beta","delta","distribution"))
    C$data$pscale <- L$data$pscale[which(parmatrix$group %in% c("alpha","gamma","beta","delta","distribution"))]
    jacfun <- MakeADFun(data = C$data, parameters = C$parameters, map = C$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    return(list(jacfun = jacfun, active_constraints_index = active_constraints_index, jacobian_matrix = jacobian_matrix))
}

.tmb_initializa_aparch <- function(spec, ...)
{
    estimate <- include <- parameter <- group <- value <- .N <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda,
                 samplen = spec$model$sample_n, initmethod = spec$model$init,
                 pscale = pscale, cmodel = cmodel,
                 model = "aparch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map))
}

.estimate_aparch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_aparch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    Q <- aparch_jac_fun(L, object$parmatrix)
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$jacfun <- Q$jacfun
    env$active_constraints_index <- Q$active_constraints_index
    env$jacobian_matrix <- Q$jacobian_matrix
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "aparch"
    env$stationarity_constraint <- stationarity_constraint
    env$parmatrix <- object$parmatrix
    env$distribution <- object$distribution
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale

    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad,
                  eval_g_ineq = aparch_ineq, eval_jac_g_ineq = aparch_ineq_jac,
                  lb = lower, ub = upper, env = env, opts = control)

    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    # we don't need special flags for egarch
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_aparch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)
    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table
    kappa_table <- scaled_solution$kappa_table
    parmatrix[parameter == "omega", value := scaled_solution$target_omega]
    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                kappa = scaled_solution$kappa,
                sigma = scaled_solution$sigma,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                kappa_summary = kappa_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    rm(tmb)
    return(out)
}

.estimate_aparch_scaled <- function(pars, H, object, control, stationarity_constraint = 0.999)
{
    value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H  <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]
    scaled_L <- .tmb_initializa_aparch(scaled_object)
    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_Q <- aparch_jac_fun(scaled_L, scaled_object$parmatrix)
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "aparch"
    scaled_env$jacfun <- scaled_Q$jacfun
    scaled_env$active_constraints_index <- scaled_Q$active_constraints_index
    scaled_env$jacobian_matrix <- scaled_Q$jacobian_matrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = aparch_ineq, eval_jac_g_ineq = aparch_ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he(scaled_sol$solution)
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) sig <- sig[-seq_len(m)]
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega
    kappa <- scaled_env$tmb$report(scaled_sol$solution)$kappa

    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]
    kappa_table <- rr["kappa",]
    rm(scaled_tmb)
    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega, kappa = kappa,
                var_initial = var_initial, sigma = sig, persistence_table = persistence_table,
                variance_target_table = variance_target_table, kappa_table = kappa_table))
}


#--------------------------------
# gjgarch model
gjrgarch_jac_fun <- function(L, parmatrix)
{
    C <- list()
    C$data$cmodel <- L$data$cmodel
    C$data$model <- "gjrgarchjacobian"
    C$parameters <- L$parameters
    C$parameters$mu <- NULL
    C$parameters$omega <- NULL
    C$parameters$xi <- NULL
    C$map <- L$map
    C$map$mu <- NULL
    C$map$omega <- NULL
    C$map$xi <- NULL
    n <- sum(parmatrix$estimate)
    jacobian_matrix <- matrix(0, ncol = sum(parmatrix$estimate), nrow = 1)
    groups <- parmatrix$group
    active_constraints_index <- which(parmatrix[estimate == 1]$group %in% c("alpha","gamma","beta","distribution"))
    C$data$pscale <- L$data$pscale[which(parmatrix$group %in% c("alpha","gamma","beta","distribution"))]
    jacfun <- MakeADFun(data = C$data, parameters = C$parameters, map = C$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    return(list(jacfun = jacfun, active_constraints_index = active_constraints_index, jacobian_matrix = jacobian_matrix))
}

.tmb_initializa_gjrgarch <- function(spec, ...)
{
    estimate <- include <- parameter <- group <- value <- .N <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda,
                 samplen = spec$model$sample_n, initmethod = spec$model$init,
                 pscale = pscale, cmodel = cmodel,
                 model = "gjrgarch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map))
}

.estimate_gjrgarch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_gjrgarch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    Q <- gjrgarch_jac_fun(L, object$parmatrix)
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$jacfun <- Q$jacfun
    env$active_constraints_index <- Q$active_constraints_index
    env$jacobian_matrix <- Q$jacobian_matrix
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "gjrgarch"
    env$stationarity_constraint <- stationarity_constraint
    env$parmatrix <- object$parmatrix
    env$distribution <- object$distribution
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale

    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad,
                  eval_g_ineq = gjrgarch_ineq, eval_jac_g_ineq = gjrgarch_ineq_jac,
                  lb = lower, ub = upper, env = env, opts = control)
    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    # we don't need special flags for gjrgarch
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_gjrgarch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)

    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table
    kappa_table <- scaled_solution$kappa_table

    parmatrix[parameter == "omega", value := scaled_solution$target_omega]
    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                kappa = scaled_solution$kappa,
                sigma = scaled_solution$sigma,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                kappa_summary = kappa_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    rm(tmb)
    return(out)
}

.estimate_gjrgarch_scaled <- function(pars, H, object, control, stationarity_constraint = 0.999)
{
    value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H  <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]
    scaled_L <- .tmb_initializa_gjrgarch(scaled_object)
    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_Q <- gjrgarch_jac_fun(scaled_L, scaled_object$parmatrix)
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "gjrgarch"
    scaled_env$jacfun <- scaled_Q$jacfun
    scaled_env$active_constraints_index <- scaled_Q$active_constraints_index
    scaled_env$jacobian_matrix <- scaled_Q$jacobian_matrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = gjrgarch_ineq, eval_jac_g_ineq = gjrgarch_ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he(scaled_sol$solution)
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) sig <- sig[-seq_len(m)]
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega
    kappa <- scaled_env$tmb$report(scaled_sol$solution)$kappa

    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]
    kappa_table <- rr["kappa",]

    rm(scaled_tmb)
    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega, kappa = kappa,
                persistence_table = persistence_table,
                variance_target_table = variance_target_table,
                kappa_table = kappa_table,
                var_initial = var_initial, sigma = sig))
}

#--------------------------------
# fgarch model
fgarch_jac_fun <- function(L, parmatrix)
{
    C <- list()
    C$data$cmodel <- L$data$cmodel
    C$data$model <- "fgarchjacobian"
    C$parameters <- L$parameters
    C$parameters$mu <- NULL
    C$parameters$omega <- NULL
    C$parameters$xi <- NULL
    C$map <- L$map
    C$map$mu <- NULL
    C$map$omega <- NULL
    C$map$xi <- NULL
    n <- sum(parmatrix$estimate)
    jacobian_matrix <- matrix(0, ncol = sum(parmatrix$estimate), nrow = 1)
    groups <- parmatrix$group
    active_constraints_index <- which(parmatrix[estimate == 1]$group %in% c("alpha","gamma","eta","beta","delta","distribution"))
    C$data$pscale <- L$data$pscale[which(parmatrix$group %in% c("alpha","gamma","eta","beta","delta","distribution"))]
    jacfun <- MakeADFun(data = C$data, parameters = C$parameters, map = C$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    return(list(jacfun = jacfun, active_constraints_index = active_constraints_index, jacobian_matrix = jacobian_matrix))
}

.tmb_initializa_fgarch <- function(spec, ...)
{
    estimate <- include <- parameter <- group <- value <- .N <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda,
                 samplen = spec$model$sample_n, initmethod = spec$model$init,
                 pscale = pscale, cmodel = cmodel,
                 model = "fgarch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map))
}

.estimate_fgarch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_fgarch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    Q <- fgarch_jac_fun(L, object$parmatrix)
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$jacfun <- Q$jacfun
    env$active_constraints_index <- Q$active_constraints_index
    env$jacobian_matrix <- Q$jacobian_matrix
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "fgarch"
    env$stationarity_constraint <- stationarity_constraint
    env$parmatrix <- object$parmatrix
    env$distribution <- object$distribution
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale

    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad,
                  eval_g_ineq = fgarch_ineq, eval_jac_g_ineq = fgarch_ineq_jac,
                  lb = lower, ub = upper, env = env, opts = control)

    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    # we don't need special flags for egarch
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_fgarch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)

    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table
    kappa_table <- scaled_solution$kappa_table

    parmatrix[parameter == "omega", value := scaled_solution$target_omega]
    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                kappa = scaled_solution$kappa,
                sigma = scaled_solution$sigma,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                kappa_summary = kappa_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    rm(tmb)
    return(out)
}

.estimate_fgarch_scaled <- function(pars, H, object, control, stationarity_constraint = 0.999)
{
    value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H  <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]
    scaled_L <- .tmb_initializa_fgarch(scaled_object)
    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_Q <- fgarch_jac_fun(scaled_L, scaled_object$parmatrix)
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "fgarch"
    scaled_env$jacfun <- scaled_Q$jacfun
    scaled_env$active_constraints_index <- scaled_Q$active_constraints_index
    scaled_env$jacobian_matrix <- scaled_Q$jacobian_matrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = fgarch_ineq, eval_jac_g_ineq = fgarch_ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he(scaled_sol$solution)
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) sig <- sig[-seq_len(m)]
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega
    kappa <- scaled_env$tmb$report(scaled_sol$solution)$kappa

    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]
    kappa_table <- rr["kappa",]


    rm(scaled_tmb)
    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega, kappa = kappa,
                persistence_table = persistence_table,
                variance_target_table = variance_target_table,
                kappa_table = kappa_table,
                var_initial = var_initial, sigma = sig))
}


#--------------------------------
# cgarch model
.tmb_initializa_cgarch <- function(spec, ...)
{
    include <- value <- group <- .N <- NULL
    parmatrix <- copy(spec$parmatrix)
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    active_constraints_index1 <- parmatrix[estimate == 1]
    active_constraints_index1 <- which(active_constraints_index1$group %in% c("rho","alpha","beta"))
    active_constraints_sign1 <- rep(1, length(active_constraints_index1))
    active_constraints_sign1[1] <- -1

    active_constraints_index2 <- parmatrix[estimate == 1]
    active_constraints_index2 <- which(active_constraints_index2$group %in% c("rho"))
    active_constraints_sign2 <- rep(1, length(active_constraints_index2))
    active_constraints_sign2[1] <- 1

    active_constraints_index3 <- parmatrix[estimate == 1]
    active_constraints_index3 <- which(active_constraints_index3$group %in% c("phi","beta"))
    active_constraints_sign3 <- rep(1, length(active_constraints_index3))
    active_constraints_sign3[-1] <- -1

    active_constraints_index4 <- parmatrix[estimate == 1]
    active_constraints_index4 <- which(active_constraints_index4$group %in% c("phi","alpha"))
    active_constraints_sign4 <- rep(1, length(active_constraints_index4))
    active_constraints_sign4[-1] <- -1

    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda, samplen = spec$model$sample_n, initmethod = spec$model$init, pscale = pscale, cmodel = cmodel, model = "cgarch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map, active_constraints_index = list(
        active_constraints_index1, active_constraints_index2,active_constraints_index3,active_constraints_index4),
        active_constraints_sign = list(active_constraints_sign1, active_constraints_sign2,active_constraints_sign3,active_constraints_sign4)))
}


.estimate_cgarch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- group <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_cgarch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, atomic = TRUE, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "cgarch"
    env$stationarity_constraint <- stationarity_constraint
    env$active_constraints_index <- L$active_constraints_index
    env$active_constraints_sign <- L$active_constraints_sign
    env$parmatrix <- object$parmatrix
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale
    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad, eval_g_ineq = cgarch_ineq, eval_jac_g_ineq = cgarch_ineq_jac, lb = lower, ub = upper, env = env, opts = control)
    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    flag <- FALSE
    if (any(pmatrix[estimate == 1 & group == "alpha"]$value == pmatrix[estimate == 1 & group == "alpha"]$lower)) {
        flag <- TRUE
        pmatrix[estimate == 1 & group == "alpha" & value == lower, value := 0.01]
        pmatrix[estimate == 1 & group == "alpha" & value == lower, lower := lower + 1e-6]
    }
    if (any(pmatrix[estimate == 1 & group == "beta"]$value == pmatrix[estimate == 1 & group == "beta"]$lower)) {
        flag <- TRUE
        pmatrix[estimate == 1 & group == "beta" & value == lower, value := 0.8]
        pmatrix[estimate == 1 & group == "beta" & value == lower, lower := lower + 1e-6]
    }
    # special case trap
    if (flag) {
        warning("\nfound potential problem in optimization with at least one parameter stuck in lower bound. Adjusting bounds and re-estimating to check.")
        lower <- pmatrix[estimate == 1]$lower * pmatrix[estimate == 1]$scale
        upper <- pmatrix[estimate == 1]$upper * 1/pmatrix[estimate == 1]$scale
        sol_check <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad, eval_g_ineq = cgarch_ineq, eval_jac_g_ineq = cgarch_ineq_jac, lb = lower, ub = upper, env = env, opts = control)
        if (sol_check$objective < sol$objective) {
            sol <- sol_check
            env$parmatrix <- pmatrix
            object$parmatrix <- pmatrix
        }
    }
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_cgarch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)

    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table

    parmatrix[parameter == "omega", value := scaled_solution$target_omega]

    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                sigma = scaled_solution$sigma,
                permanent_component = scaled_solution$permanent_component,
                transitory_component = scaled_solution$sigma - scaled_solution$permanent_component,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    return(out)
}

.estimate_cgarch_scaled <- function(pars, H, object, control, stationarity_constraint)
{
    parameter <- value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]

    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_L <- .tmb_initializa_cgarch(scaled_object)
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "cgarch"
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$active_constraints_index <- scaled_L$active_constraints_index
    scaled_env$active_constraints_sign <- scaled_L$active_constraints_sign
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_ineq = cgarch_ineq, eval_jac_g_ineq = cgarch_ineq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he()
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    permanent_component <- scaled_env$tmb$report(scaled_sol$solution)$permanent_component
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) {
        sig <- sig[-seq_len(m)]
        permanent_component <- permanent_component[-seq_len(m)]
    }
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega

    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]


    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega, permanent_component = permanent_component,
                persistence_table = persistence_table, variance_target_table = variance_target_table,
                var_initial = var_initial, sigma = sig))
}



#--------------------------------
# igarch model
.tmb_initializa_igarch <- function(spec, ...)
{
    include <- value <- group <- .N <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    pscale <- parmatrix$scale
    active_constraints_index <- parmatrix[estimate == 1]
    active_constraints_index <- which(active_constraints_index$group %in% c("alpha","beta"))
    data <- list(y = y, v = v, backcast_lambda = spec$model$backcast_lambda, samplen = spec$model$sample_n, initmethod = spec$model$init, pscale = pscale, cmodel = cmodel, model = "garch")
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map, active_constraints_index = active_constraints_index))
}


.estimate_igarch <- function(object, solver, control, stationarity_constraint = 0.999, ...)
{
    parameter <- group <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    L <- .tmb_initializa_igarch(object)
    tmb <- MakeADFun(data = L$data, parameters = L$parameters, atomic = TRUE, map = L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    env <- new.env()
    env$fun <- L$fun
    env$grad <- L$grad
    env$hess <- L$hess
    env$tmb <- tmb
    env$pindex <- L$pindex
    env$llh <- 1
    env$model <- "garch"
    env$stationarity_constraint <- stationarity_constraint
    env$active_constraints_index <- L$active_constraints_index
    env$parmatrix <- object$parmatrix
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale
    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad, eval_g_eq = igarch_eq, eval_jac_g_eq = igarch_eq_jac,
                  lb = lower, ub = upper, env = env, opts = control)
    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]
    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_igarch_scaled(optimal_pars, H, object, control, stationarity_constraint = stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution, scaled_solution$env$fun, scaled_solution$env$grad, scaled_solution$env$hess, scaled_solution$lower, scaled_solution$upper, scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)

    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table

    parmatrix[parameter == "omega", value := scaled_solution$target_omega]
    out <- list(parmatrix = parmatrix, scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                sigma = scaled_solution$sigma,
                loglik = scaled_solution$solution$objective,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec)
    class(out) <- "tsgarch.estimate"
    return(out)
}

.estimate_igarch_scaled <- function(pars, H, object, control, stationarity_constraint)
{
    parameter <- value <- NULL
    test <- .is_positive_definite(H)
    if (!test) H <- .make_positive_definite(H)
    par_scale <- sqrt(1/diag(H))
    # stage 2: scaled estimation
    scaled_object <- object
    # make a copy of parmatrix rather than leave the pointer in place
    scaled_object$parmatrix <- copy(object$parmatrix)
    scaled_object$parmatrix[estimate == 1, scale := par_scale]
    scaled_object$parmatrix[estimate == 1, value := pars / par_scale]

    scaled_init_pars <- pars *  1/par_scale
    scaled_lower <- scaled_object$parmatrix[estimate == 1]$lower * 1/par_scale
    scaled_upper <-  scaled_object$parmatrix[estimate == 1]$upper * 1/par_scale
    scaled_L <- .tmb_initializa_igarch(scaled_object)
    scaled_tmb <- MakeADFun(data = scaled_L$data, parameters = scaled_L$parameters, map = scaled_L$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_env <- new.env()
    scaled_env$fun <- scaled_L$fun
    scaled_env$grad <- scaled_L$grad
    scaled_env$hess <- scaled_L$hess
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- "garch"
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_env$active_constraints_index <- scaled_L$active_constraints_index
    scaled_sol <- nloptr(x0 = scaled_init_pars, eval_f = scaled_env$fun, eval_grad_f = scaled_env$grad,
                         eval_g_eq = igarch_eq, eval_jac_g_eq = igarch_eq_jac,
                         lb = scaled_lower, ub = scaled_upper, env = scaled_env, opts = control)
    scaled_sol$par_scale <- par_scale
    hessian <- scaled_tmb$he()
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    if (m > 0) {
        sig <- sig[-seq_len(m)]
    }
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega

    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]


    return(list(solution = scaled_sol, env = scaled_env, lower = scaled_lower, upper = scaled_upper,
                hessian = hessian, scores = scores, par_scale = par_scale,
                target_omega = target_omega, persistence_table = persistence_table,
                variance_target_table = variance_target_table,
                var_initial = var_initial, sigma = sig))
}
