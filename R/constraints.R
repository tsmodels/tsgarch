# nloptr constraints on persistence < 1 (we use 0.999 since g(x) <= 0)
garch_ineq <- function(pars, env)
{
    p <- .persistence(pars, env)
    # g(x) <= 0
    return(p - env$stationarity_constraint)
}

# garch constraint Jacobian
garch_ineq_jac <- function(pars, env)
{
    parmatrix <- env$parmatrix
    p <- parmatrix[estimate == 1]$value * 0
    p[env$active_constraints_index] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index]
    return(p)
}

# egarch constraint Jacobian
egarch_ineq_jac <- function(pars, env)
{
    parmatrix <- env$parmatrix
    p <- parmatrix[estimate == 1]$value * 0
    p[env$active_constraints_index] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index]
    return(p)
}

# aparch ineq jacobian --------------------------------------------------------
aparch_ineq <- function(pars, env)
{
    p <- env$jacfun$fn(pars[env$active_constraints_index])
    return(p - env$stationarity_constraint)
}

aparch_ineq_jac <- function(pars, env)
{
    out <- env$jacobian_matrix * 0
    jacs <- env$jacfun$gr(pars[env$active_constraints_index])
    out[,env$active_constraints_index] <- jacs
    return(out)
}

# gjrgarch ineq jacobian --------------------------------------------------------
gjrgarch_ineq <- function(pars, env)
{
    p <- env$jacfun$fn(pars[env$active_constraints_index])
    return(p - env$stationarity_constraint)
}

gjrgarch_ineq_jac <- function(pars, env)
{
    out <- env$jacobian_matrix * 0
    jacs <- env$jacfun$gr(pars[env$active_constraints_index])
    out[,env$active_constraints_index] <- jacs
    return(out)
}

# fgarch ineq jacobian --------------------------------------------------------
fgarch_ineq <- function(pars, env)
{
    p <- env$jacfun$fn(pars[env$active_constraints_index])
    return(p - env$stationarity_constraint)
}

fgarch_ineq_jac <- function(pars, env)
{
    out <- env$jacobian_matrix * 0
    jacs <- env$jacfun$gr(pars[env$active_constraints_index])
    out[,env$active_constraints_index] <- jacs
    return(out)
}


# cgarch ineq jacobian --------------------------------------------------------
cgarch_ineq <- function(pars, env)
{
    # return the transitory component persistence
    group <- value <- NULL
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
    # phi < alpha
    p4 <-  phi - sum(parmatrix[group == "alpha"]$value * parmatrix[group == "alpha"]$scale)
    return(c(p1, p2, p3, p4))
}

cgarch_ineq_jac <- function(pars, env)
{
    parmatrix <- env$parmatrix
    p <- rbind(parmatrix[estimate == 1]$value * 0, parmatrix[estimate == 1]$value * 0, parmatrix[estimate == 1]$value * 0, parmatrix[estimate == 1]$value * 0)
    p[1, env$active_constraints_index[[1]]] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index[[1]]] * env$active_constraints_sign[[1]]
    p[2, env$active_constraints_index[[2]]] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index[[2]]] * env$active_constraints_sign[[2]]
    p[3, env$active_constraints_index[[3]]] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index[[3]]] * env$active_constraints_sign[[3]]
    p[4, env$active_constraints_index[[4]]] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index[[4]]] * env$active_constraints_sign[[4]]
    return(p)
}

# igarch eq jacobian --------------------------------------------------------
igarch_eq <- function(pars, env)
{
    p <- .persistence(pars, env)
    # g(x) == 0
    return(p - 1.0)
}

# igarch constraint Jacobian
igarch_eq_jac <- function(pars, env)
{
    parmatrix <- env$parmatrix
    p <- parmatrix[estimate == 1]$value * 0
    p[env$active_constraints_index] <- env$parmatrix[estimate == 1]$scale[env$active_constraints_index]
    return(p)
}
