# !diagnostics suppress=copy

#--------------------------------
# common wrappers for the objective, gradient and hessian functions
garch_fun <- function(pars, env)
{
    estimate <- parameter <- value <- NULL
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

solve_model <- function(init_pars, env, const, lower, upper, control) {
    sol <- nloptr(x0 = init_pars, eval_f = env$fun, eval_grad_f = env$grad,
                  eval_g_ineq = const$inequality_constraint,
                  eval_jac_g_ineq = const$inequality_jacobian,
                  eval_g_eq = const$equality_constraint,
                  eval_jac_g_eq = const$equality_jacobian,
                  lb = lower, ub = upper, env = env, opts = control)
    return(sol)
}


# solve_solnp_model <- function(init_pars, env, const, lower, upper, control) {
#     m <- length(const$inequality_constraint(init_pars, env))
#     ineq_lower <- rep(-1000, m)
#     ineq_upper <- rep(0, m)
#     sol <- Rsolnp::csolnp(pars = init_pars, fn = env$fun, gr = env$grad,
#                   ineq_fn = const$inequality_constraint,
#                   ineq_lower = ineq_lower, ineq_upper = ineq_upper,
#                   ineq_jac = const$inequality_jacobian,
#                   eq_fn = const$equality_constraint,
#                   eq_jac = const$equality_jacobian,
#                   lower = lower, upper = upper, env = env, control = list(trace = 1))
#     sol$solution <- sol$pars
#     return(sol)
# }

# common functions


.tmb_initialize_model <- function(spec, ...)
{
    estimate <- include <- value <- group <- .N <- NULL
    # copy rather than alias: this function mutates parmatrix (both the
    # pre-existing "include" helper column below and, when an entire ar/ma
    # polynomial is fixed, the arpacf/mapacf values themselves - see below)
    # and must not leak those changes back into the caller's spec object.
    parmatrix <- copy(spec$parmatrix)
    # A spec serialized by a version predating the ARMA-X mean equation has no
    # arpacf/mapacf/tau rows at all. MakeADFun would otherwise fail on the
    # first missing parameter with an opaque "Error when reading the variable:
    # 'arpacf'", so report what is actually wrong and how to recover.
    if (!all(c("arpacf","mapacf","tau") %in% parmatrix$group)) {
        stop("\nspec$parmatrix has no mean equation (arpacf/mapacf/tau) rows: this specification was created by an older version of tsgarch. Re-create it with garch_modelspec().", call. = FALSE)
    }
    arma_order <- spec$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    # Users can fix an entire AR and/or MA polynomial at target coefficients
    # simply by setting value = the desired ar/ma coefficients (not the
    # internal pacf parameterization) and estimate = 0 on every row of the
    # arpacf/mapacf group (see arma_coefficients()/garch_modelspec()); this
    # is the one place that "magic" needs to be resolved into the actual
    # raw pacf parameters TMB expects (which will, internally, re-apply the
    # forward pacf_to_ar/pacf_to_ma transform and recover exactly the fixed
    # target, since ar_to_pacf/ma_to_pacf are its exact inverse).
    # arma_block_status() rejects partial (mixed estimate) fixing within a
    # group with an informative error.
    if (arma_order[1] > 0 && identical(arma_block_status(parmatrix, "arpacf"), "fixed")) {
        target_ar <- parmatrix[group == "arpacf"]$value
        validate_arma_fixed_target(target_ar, "ar")
        parmatrix[group == "arpacf", value := ar_to_pacf(target_ar)]
    }
    if (arma_order[2] > 0 && identical(arma_block_status(parmatrix, "mapacf"), "fixed")) {
        target_ma <- parmatrix[group == "mapacf"]$value
        validate_arma_fixed_target(target_ma, "ma")
        parmatrix[group == "mapacf", value := ma_to_pacf(target_ma)]
    }
    parmatrix[,include := 1]
    parmatrix[estimate == 0, include := as.numeric(NA)]
    map <- lapply(split(parmatrix[,list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    parameters <- lapply(split(parmatrix[,list(value, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$value))
    cmodel <- spec$model_options
    # Specs serialized by an older version carry a shorter model_options. Every
    # released version up to 1.0.4 wrote six elements ([maxpq, arch, garch,
    # variance_targeting, multiplicative, distribution]) and the ar/ma/armax
    # flags were appended only afterwards, so the shortfall is not always one
    # element. The templates read cmodel(6..8) unconditionally, so pad to the
    # full length rather than by a fixed count.
    if (length(cmodel) < 9L) cmodel <- c(cmodel, rep(0L, 9L - length(cmodel)))
    # augment data with max(p,q) vectors
    y <- c(rep(0, cmodel[1]), as.numeric(spec$target$y))
    v <- spec$vreg$vreg
    v <- rbind(matrix(0, ncol = ncol(v), nrow = cmodel[1]), v)
    x <- spec$xreg$xreg
    if (is.null(x)) x <- matrix(0, ncol = 1, nrow = NROW(spec$target$y))
    x <- rbind(matrix(0, ncol = ncol(x), nrow = cmodel[1]), x)
    pscale <- parmatrix$scale
    model <- spec$model$model
    if (model == "igarch") {
        model <- "garch"
    }
    data <- list(y = y, v = v, x = x, backcast_lambda = spec$model$backcast_lambda, samplen = spec$model$sample_n, initmethod = spec$model$init, pscale = pscale, cmodel = cmodel, model = model)
    return(list(fun = garch_fun, grad = garch_grad, hess = garch_hess, data = data, parameters = parameters, map = map))
}

.estimate_garch_model <- function(object, solver, control, stationarity_constraint = 0.999, keep_tmb = FALSE, ...)
{
    estimate <- parameter <- value <- NULL
    solver <- match.arg(solver, choices = "nloptr")
    model_init <- .tmb_initialize_model(object)
    const <- setup_garch_constraints(object, model_init)
    tmb <- MakeADFun(data = model_init$data, parameters = model_init$parameters, map = model_init$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    env <- new.env()
    env$fun <- model_init$fun
    env$grad <- model_init$grad
    env$hess <- model_init$hess
    env$tmb <- tmb
    env$llh <- 1
    env$model <- object$model$model
    env$stationarity_constraint <- stationarity_constraint
    env$active_constraints <- const$active_constraints
    env$jacobian_matrix <- const$jacobian_matrix
    env$jacfun <- const$jacfun
    env$parmatrix <- object$parmatrix
    env$distribution <- object$distribution
    init_pars <- object$parmatrix[estimate == 1]$value * 1/object$parmatrix[estimate == 1]$scale
    lower <- object$parmatrix[estimate == 1]$lower * 1/object$parmatrix[estimate == 1]$scale
    upper <- object$parmatrix[estimate == 1]$upper * 1/object$parmatrix[estimate == 1]$scale
    sol <- solve_model(init_pars = init_pars, env = env, const = const, lower = lower, upper = upper, control = control)
    #sol <- solve_solnp_model(init_pars = init_pars, env = env, const = const, lower = lower, upper = upper, control = control)
    pmatrix <- copy(env$parmatrix)
    pmatrix[estimate == 1, value := sol$solution]

    optimal_pars <- sol$solution
    # check hessian and use for scaling
    H <- tmb$he(optimal_pars)
    if (any(is.na(H))) {
        warning("\nunable to calculate hessian for parameter scaling in initial step. Reverting to numerical estimation.")
        H <- hessian(env$fun, x = optimal_pars, env = env)
    }
    object$parmatrix <- pmatrix
    scaled_solution <- .estimate_garch_model_scaled(optimal_pars, H, object, control, stationarity_constraint)
    D <- solve(diag(scaled_solution$par_scale, length(scaled_solution$par_scale), length(scaled_solution$par_scale)))
    pars <- scaled_solution$solution$solution * scaled_solution$solution$par_scale
    solve_conditions <- solver_conditions(scaled_solution$solution$solution,
                                          scaled_solution$env$fun,
                                          scaled_solution$env$grad,
                                          scaled_solution$env$hess,
                                          scaled_solution$lower,
                                          scaled_solution$upper,
                                          scaled_solution$env)
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    spec <- object
    spec$parmatrix <- NULL
    spec$model$var_initial <- scaled_solution$var_initial
    conditional_mu <- scaled_solution$conditional_mu
    # the variance target is the unconditional variance of the ARMA
    # innovations eps = y - conditional_mu over the full sample (see
    # vignettes/garch_models.Rmd, "Variance Targeting"); fall back to the
    # constant-mean deviation when no ARMA mean equation is present
    if (!is.null(conditional_mu)) {
        constant_variance <- mean((as.numeric(object$target$y_orig) - conditional_mu)^2)
    } else {
        constant_variance <- mean((object$target$y_orig - (parmatrix[parameter == "mu"]$value * parmatrix[parameter == "mu"]$scale))^2)
    }
    persistence_table <- scaled_solution$persistence_table
    variance_target_table <- scaled_solution$variance_target_table
    kappa_table <- scaled_solution$kappa_table
    kappa <- scaled_solution$kappa
    permanent_component = scaled_solution$permanent_component
    transitory_component = scaled_solution$transitory_component

    parmatrix[parameter == "omega", value := scaled_solution$target_omega]

    out <- list(parmatrix = parmatrix,
                scaled_hessian = scaled_solution$hessian,
                scaled_scores = scaled_solution$scores,
                parameter_scale = scaled_solution$par_scale,
                conditions = solve_conditions,
                var_initial = scaled_solution$var_initial,
                arch_initial = scaled_solution$arch_initial,
                constant_variance = constant_variance,
                target_omega = scaled_solution$target_omega,
                sigma = scaled_solution$sigma,
                permanent_component = permanent_component,
                transitory_component = transitory_component,
                loglik = scaled_solution$solution$objective,
                lik_vector = scaled_solution$ll_vector,
                nobs = NROW(spec$target$y),
                persistence_summary = persistence_table,
                variance_target_summary = variance_target_table,
                kappa_summary = kappa_table,
                kappa = kappa,
                # extra degree of freedom for the init_variance
                npars = NROW(parmatrix[estimate == 1]) + 1,
                spec = spec,
                conditional_mu = conditional_mu,
                arma_summary = scaled_solution$arma_table)
    if (keep_tmb) out$tmb <- tmb
    class(out) <- "tsgarch.estimate"
    return(out)
}

.estimate_garch_model_scaled <- function(pars, H, object, control, stationarity_constraint)
{
    estimate <- parameter <- value <- NULL
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
    scaled_model_init <- .tmb_initialize_model(scaled_object)
    scaled_const <- setup_garch_constraints(scaled_object, scaled_model_init)
    scaled_tmb <- MakeADFun(data = scaled_model_init$data, parameters = scaled_model_init$parameters, map = scaled_model_init$map, silent = TRUE, DLL = "tsgarch_TMBExports")
    scaled_env <- new.env()
    scaled_env$fun <- scaled_model_init$fun
    scaled_env$grad <- scaled_model_init$grad
    scaled_env$hess <- scaled_model_init$hess
    scaled_env$active_constraints <- scaled_const$active_constraints
    scaled_env$jacobian_matrix <- scaled_const$jacobian_matrix
    scaled_env$jacfun <- scaled_const$jacfun
    scaled_env$tmb <- scaled_tmb
    scaled_env$llh <- 1
    scaled_env$model <- object$model$model
    scaled_env$distribution <- object$distribution
    scaled_env$parmatrix <- scaled_object$parmatrix
    scaled_env$stationarity_constraint <- stationarity_constraint
    scaled_sol <- solve_model(init_pars = scaled_init_pars, env = scaled_env, const = scaled_const, lower = scaled_lower, upper = scaled_upper, control = control)
    #scaled_sol <- solve_solnp_model(init_pars = scaled_init_pars, env = scaled_env, const = scaled_const, lower = scaled_lower, upper = scaled_upper, control = control)
    scaled_sol$par_scale <- par_scale
    # evaluate at the solution actually returned rather than relying on TMB's
    # default (x = last.par, the optimizer's last evaluated point, which need
    # not be the point it returns); this hessian feeds bread() and therefore
    # every standard error
    hessian <- scaled_tmb$he(scaled_sol$solution)
    if (any(is.na(hessian))) {
        warning("\nunable to calculate hessian for parameter scaling in scaling step. Reverting to numerical estimation.")
        hessian <- hessian(scaled_env$fun, x = scaled_sol$solution, env = scaled_env)
    }
    scores <- jacobian(score_function, scaled_sol$solution, env = scaled_env)
    m <- object$model_options[1]
    sig <- scaled_env$tmb$report(scaled_sol$solution)$sigma
    var_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_variance
    arch_initial <- scaled_env$tmb$report(scaled_sol$solution)$initial_arch
    target_omega <- scaled_env$tmb$report(scaled_sol$solution)$target_omega
    rr <- summary(sdreport(scaled_tmb, par.fixed = scaled_sol$solution, getReportCovariance = T), p.value = TRUE)
    persistence_table <- rr["persistence", ]
    variance_target_table <- rr["target_omega", ]
    ll_vector <- -1 * log(scaled_env$tmb$report(scaled_sol$solution)$ll_vector)
    if (object$model$model %in% c("gjrgarch","egarch","aparch","fgarch")) {
        kappa_table <- scaled_env$tmb$report(scaled_sol$solution)$kappa_table
        kappa_table <- rr["kappa",]
        kappa <- scaled_env$tmb$report(scaled_sol$solution)$kappa
    } else {
        kappa_table <- NULL
        kappa <- 1
    }
    if (object$model$model %in% c("cgarch")) {
        permanent_component <- scaled_env$tmb$report(scaled_sol$solution)$permanent_component
        transitory_component <- scaled_env$tmb$report(scaled_sol$solution)$transitory_component
        if (m > 0) {
            transitory_component <- transitory_component[-seq_len(m)]
            permanent_component <- permanent_component[-seq_len(m)]
        }
    } else{
        permanent_component <- NULL
        transitory_component <- NULL
    }
    if (m > 0) {
        sig <- sig[-seq_len(m)]
    }
    # conditional_mean(t) is the model's fitted conditional mean of y at
    # time t (mu + the AR/MA deviation term), reported by the "garch" TMB
    # template only when arma > c(0,0) (see garchfun.hpp); every other case
    # falls back to NULL so that fitted.tsgarch.estimate() keeps using its
    # prior rep(mu, n) behavior. residuals are always y - conditional_mu
    # (see residuals.tsgarch.estimate()), so no separate residual field is
    # needed here.
    conditional_mu <- NULL
    arma_table <- NULL
    arma_order <- object$model$arma
    if (is.null(arma_order)) arma_order <- c(0,0)
    # also capture the conditional mean when the model has mean regressors
    # with no ARMA dynamics: conditional_mean(i) = mu + xtau(i) is then
    # time-varying and fitted()/residuals() need it
    if (sum(arma_order) > 0 || isTRUE(object$xreg$include_xreg)) {
        # sig has already been trimmed to the actual (non-padded) length above
        conditional_mu <- tail(scaled_env$tmb$report(scaled_sol$solution)$conditional_mean, length(sig))
    }
    if (sum(arma_order) > 0) {
        ar_order <- object$model$arma[1]
        ma_order <- object$model$arma[2]
        arma_table <- list(
            ar = if (ar_order > 0) rr[rownames(rr) == "arma_ar", , drop = FALSE][seq_len(ar_order), , drop = FALSE] else NULL,
            ma = if (ma_order > 0) rr[rownames(rr) == "arma_ma", , drop = FALSE][seq_len(ma_order), , drop = FALSE] else NULL)
    }
    rm(scaled_tmb)
    return(list(solution = scaled_sol,
                env = scaled_env,
                lower = scaled_lower,
                upper = scaled_upper,
                hessian = hessian,
                scores = scores,
                par_scale = par_scale,
                target_omega = target_omega,
                persistence_table = persistence_table,
                variance_target_table = variance_target_table,
                var_initial = var_initial,
                arch_initial = arch_initial,
                kappa_table = kappa_table,
                kappa = kappa,
                sigma = sig,
                permanent_component = permanent_component,
                transitory_component = transitory_component,
                ll_vector = ll_vector,
                conditional_mu = conditional_mu,
                arma_table = arma_table
                ))
}
