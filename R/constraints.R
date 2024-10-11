# !diagnostics suppress=copy

# nloptr constraints on persistence < 1 (we use 0.999 since g(x) <= 0)
# g(x) <= 0

active_constraints_group <- function(model)
{
    # if any of these are active, then the constraint or part of it is also active
    switch(model,
           "garch" = c("alpha", "beta"),
           "igarch" = c("alpha","beta"),
           "egarch" = c("beta"),
           "gjrgarch" = c("alpha", "gamma", "beta", "distribution"),
           "aparch" = c("alpha", "beta", "gamma", "delta", "distribution"),
           "fgarch" = c("alpha", "beta", "gamma", "eta", "delta", "distribution"),
           "cgarch" = c("alpha","beta","phi","rho")
    )
}

setup_garch_constraints <- function(spec, model_init, ...)
{
    estimate <- NULL
    model <- spec$model$model
    .group <- active_constraints_group(model)
    jacfun <- NULL
    jacobian_matrix <- NULL
    inequality_constraint <- NULL
    inequality_jacobian <- NULL
    equality_constraint <- NULL
    equality_jacobian <- NULL
    active_constraints <- NULL
    estimeable_parameters <- spec$parmatrix[estimate == 1]
    n_pars <- NROW(estimeable_parameters)
    persistence_constraint <- which(estimeable_parameters$group %in% .group)
    if (length(persistence_constraint) > 0) {
        if (model %in% c("garch","egarch")) {
            active_constraints <- persistence_constraint
            inequality_constraint <- linear_inequality
            inequality_jacobian <- linear_inequality_jacobian
            jacobian_matrix <- matrix(0, ncol = n_pars, nrow = 1)
        } else if (model %in% c("aparch","fgarch")) {
            inequality_constraint <- nonlinear_inequality
            inequality_jacobian <- nonlinear_inequality_jacobian
            tmb_jac <- list()
            tmb_jac$data$cmodel <- model_init$data$cmodel
            if (model == "aparch") {
                tmb_jac$data$model <- "aparchjacobian"
            } else if (model == "fgarch") {
                tmb_jac$data$model <- "fgarchjacobian"
            }
            tmb_jac$parameters <- model_init$parameters
            tmb_jac$parameters$mu <- NULL
            tmb_jac$parameters$omega <- NULL
            tmb_jac$parameters$xi <- NULL
            tmb_jac$map <- model_init$map
            tmb_jac$map$mu <- NULL
            tmb_jac$map$omega <- NULL
            tmb_jac$map$xi <- NULL
            tmb_jac$data$pscale <- model_init$data$pscale[which(spec$parmatrix$group %in% .group)]
            jacfun <- MakeADFun(data = tmb_jac$data, parameters = tmb_jac$parameters, map = tmb_jac$map, silent = TRUE, DLL = "tsgarch_TMBExports")
            jacobian_matrix <- matrix(0, ncol = n_pars, nrow = 1)
            active_constraints <- persistence_constraint
        } else if (model == "gjrgarch") {
            # multiple constraints
            inequality_constraint <- gjr_inequality
            inequality_jacobian <- gjr_inequality_jacobian
            tmb_jac <- list()
            tmb_jac$data$cmodel <- model_init$data$cmodel
            tmb_jac$data$model <- "gjrgarchjacobian"
            tmb_jac$parameters <- model_init$parameters
            tmb_jac$parameters$mu <- NULL
            tmb_jac$parameters$omega <- NULL
            tmb_jac$parameters$xi <- NULL
            tmb_jac$map <- model_init$map
            tmb_jac$map$mu <- NULL
            tmb_jac$map$omega <- NULL
            tmb_jac$map$xi <- NULL
            tmb_jac$data$pscale <- model_init$data$pscale[which(spec$parmatrix$group %in% .group)]
            jacfun <- MakeADFun(data = tmb_jac$data, parameters = tmb_jac$parameters, map = tmb_jac$map, silent = TRUE, DLL = "tsgarch_TMBExports")
            active_constraints <- vector(mode = "list", length = 1 + model_init$data$cmodel[2])
            active_constraints[[1]] <- persistence_constraint
            jacobian_matrix <- matrix(0, ncol = n_pars, nrow = 1 + model_init$data$cmodel[2])
            # pairwise constraints \alpha_j + \gamma_j > 0
            for (j in 1:model_init$data$cmodel[2]) {
                bounds_constraints <- which(estimeable_parameters$parameter %in% c(paste0("alpha",j),paste0("gamma",j)))
                if (length(bounds_constraints) > 0) {
                    active_constraints[[j + 1]] <- bounds_constraints
                    jacobian_matrix[j + 1, bounds_constraints] <- -1
                }
            }
        } else if (model == "cgarch") {
            inequality_constraint <- cgarch_inequality
            inequality_jacobian <- cgarch_inequality_jacobian
            active_constraints <- vector(mode = "list", length = 3)
            jacobian_matrix <- matrix(0, ncol = n_pars, nrow = 3)
            constraint_1 <- which(estimeable_parameters$group %in% c("rho","alpha","beta"))
            j <- 1
            if (length(constraint_1) > 0) {
                active_constraints[[j]] <- constraint_1
                rho_include <- which(estimeable_parameters$group %in% c("rho"))
                if (length(rho_include) > 0) {
                    jacobian_matrix[j,rho_include] <- -1
                }
                alpha_include <- which(estimeable_parameters$group %in% c("alpha"))
                if (length(alpha_include) > 0) {
                    jacobian_matrix[j,alpha_include] <- 1
                }
                beta_include <- which(estimeable_parameters$group %in% c("beta"))
                if (length(beta_include) > 0) {
                    jacobian_matrix[j,beta_include] <- 1
                }
            }
            j <- 2
            constraint_2 <- which(estimeable_parameters$group %in% c("rho"))
            if (length(constraint_2) > 0) {
                active_constraints[[j]] <- constraint_2
                jacobian_matrix[j, constraint_2] <- 1
            }
            j <- 3
            constraint_3 <- which(estimeable_parameters$group %in% c("phi","beta"))
            if (length(constraint_3) > 0) {
                active_constraints[[j]] <- constraint_3
                phi_include <- which(estimeable_parameters$group %in% c("phi"))
                if (length(phi_include) > 0) {
                    jacobian_matrix[j,phi_include] <- 1
                }
                beta_include <- which(estimeable_parameters$group %in% c("beta"))
                if (length(beta_include) > 0) {
                    jacobian_matrix[j,beta_include] <- -1
                }
            }
        } else if (model == "igarch") {
            active_constraints <- persistence_constraint
            equality_constraint <- linear_equality
            equality_jacobian <- linear_equality_jacobian
            jacobian_matrix <- matrix(0, ncol = n_pars, nrow = 1)
        }
    } else {
        jacobian_matrix <- matrix(0, ncol = n_pars, nrow = 1)
    }
    return(list(jacfun = jacfun, jacobian_matrix = jacobian_matrix, active_constraints = active_constraints,
                inequality_constraint = inequality_constraint, inequality_jacobian = inequality_jacobian,
                equality_constraint = equality_constraint, equality_jacobian = equality_jacobian))
}


linear_inequality <- function(pars, env)
{
    p <- .persistence(pars, env)
    return(p - env$stationarity_constraint)
}

linear_inequality_jacobian <- function(pars, env)
{
    estimate <- NULL
    jmat <- env$jacobian_matrix * 0
    jmat[1,env$active_constraints] <- env$parmatrix[estimate == 1]$scale[env$active_constraints]
    return(jmat)
}

nonlinear_inequality <- function(pars, env)
{
    # persistence calculated using TMB to obtain the Jacobian
    p <- env$jacfun$fn(pars[env$active_constraints])
    return(p - env$stationarity_constraint)
}


nonlinear_inequality_jacobian <- function(pars, env)
{
    jmat <- env$jacobian_matrix * 0
    jac <- env$jacfun$gr(pars[env$active_constraints])
    jmat[1,env$active_constraints] <- jac
    return(jmat)
}

gjr_pairwise_inequality <- function(pars, env)
{
    value <- estimate <- group <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    parmatrix$value <- parmatrix$value * parmatrix$scale
    pairwise_sums <- parmatrix[group == "alpha"]$value + parmatrix[group == "gamma"]$value
    return(pairwise_sums)
}

gjr_inequality <- function(pars, env)
{
    # the GJR model has a nonlinear persistence inequality as well as a pairwise
    # linear inequality for the pairwise sum of alpha_j and gamma_j > 0 since
    # gamma_j can be negative
    p1 <- env$jacfun$fn(pars[env$active_constraints[[1]]])
    p2 <- gjr_pairwise_inequality(pars, env)
    return(c(p1 - env$stationarity_constraint, 0 - p2))
}

gjr_inequality_jacobian <- function(pars, env)
{
    estimate <- NULL
    jmat <- env$jacobian_matrix
    jac <- env$jacfun$gr(pars[env$active_constraints[[1]]])
    jmat[1,env$active_constraints[[1]]] <- jac
    n <- length(env$active_constraints) - 1
    for (j in seq_len(n)) {
        if (length(env$active_constraints[[j]]) > 0) {
            # jmat holds the correct signs for each constraint
            jmat[j + 1,env$active_constraints[[j + 1]]] <- jmat[j + 1,env$active_constraints[[j + 1]]] * env$parmatrix[estimate == 1]$scale[env$active_constraints[[j + 1]]]
        }
    }
    return(jmat)
}

cgarch_inequality <- function(pars, env)
{
    # return the transitory component persistence
    estimate <- group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    # alpha + beta < rho
    p <- .persistence(pars, env)
    p1 <- p[2] - p[1]
    # rho < 1
    p2 <- p[1] - env$stationarity_constraint
    # phi < beta
    phi <- (parmatrix[group == "phi"]$value * parmatrix[group == "phi"]$scale)
    p3 <-  phi - sum(parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale)
    return(c(p1, p2, p3))
}

cgarch_inequality_jacobian <- function(pars, env)
{
    estimate <- NULL
    parmatrix <- env$parmatrix
    active_constraints <- env$active_constraints
    jmat <- env$jacobian_matrix
    for (j in 1:3) {
        if (length(active_constraints[[j]]) > 0) {
            jmat[j,active_constraints[[j]]] <- jmat[j,active_constraints[[j]]] * parmatrix[estimate == 1]$scale[active_constraints[[j]]]
        }
    }
    return(jmat)
}

# linear equality (igarch)
# nloptr:  g(x) == 0
linear_equality <- function(pars, env)
{
    p <- .persistence(pars, env)
    return(p - 1.0)
}

linear_equality_jacobian <- function(pars, env)
{
    estimate <- NULL
    parmatrix <- env$parmatrix
    p <- parmatrix[estimate == 1]$value * 0
    p[env$active_constraints] <- env$parmatrix[estimate == 1]$scale[env$active_constraints]
    return(p)
}
