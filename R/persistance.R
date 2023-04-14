.persistence <- function(pars, env)
{
    p <- switch(env$model,
           "garch" = .persistence_garch(pars, env),
           "egarch" = .persistence_egarch(pars, env),
           "gjrgarch" = .persistence_gjrgarch(pars, env),
           "aparch" = .persistence_aparch(pars, env),
           "fgarch" = .persistence_fgarch(pars, env),
           "cgarch" = .persistence_cgarch(pars, env),
           # technically 1, but we calculate it for validation
           "igarch" = .persistence_garch(pars, env)
    )
    return(p)
}


.persistence_garch <- function(pars, env)
{
    parameter <- group <- value <- NULL
    parmatrix <- env$parmatrix
    parmatrix[estimate == 1, value := pars]
    p <- sum(parmatrix[group == "alpha"]$value * parmatrix[group == "alpha"]$scale) + sum(parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale)
    return(p)
}

.persistence_egarch <- function(pars, env)
{
    parameter <- group <- value <- NULL
    parmatrix <- env$parmatrix
    parmatrix[estimate == 1, value := pars]
    p <- sum(parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale)
    return(p)
}

.persistence_aparch <- function(pars, env)
{
    parameter <- group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    beta <- sum(parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale)
    skew <- parmatrix[parameter == "skew"]$value * parmatrix[parameter == "skew"]$scale
    shape <- parmatrix[parameter == "shape"]$value * parmatrix[parameter == "shape"]$scale
    lambda <- parmatrix[parameter == "lambda"]$value * parmatrix[parameter == "lambda"]$scale
    delta <- parmatrix[parameter == "delta"]$value * parmatrix[parameter == "delta"]$scale
    alpha <- parmatrix[group == "alpha"]$value * parmatrix[group == "alpha"]$scale
    gamma <- parmatrix[group == "gamma"]$value * parmatrix[group == "gamma"]$scale
    kappa <- aparch_moment_v(distribution = env$distribution, gamma = gamma, delta = delta,
                              skew = skew, shape = shape, lambda = lambda)
    alpha <- sum(alpha * kappa)
    p <- beta + alpha
    return(p)
}

.persistence_fgarch <- function(pars, env)
{
    parameter <- group <- value <- NULL
    parmatrix <- copy(env$parmatrix)
    parmatrix[estimate == 1, value := pars]
    beta <- sum(parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale)
    skew <- parmatrix[parameter == "skew"]$value * parmatrix[parameter == "skew"]$scale
    shape <- parmatrix[parameter == "shape"]$value * parmatrix[parameter == "shape"]$scale
    lambda <- parmatrix[parameter == "lambda"]$value * parmatrix[parameter == "lambda"]$scale
    delta <- parmatrix[parameter == "delta"]$value * parmatrix[parameter == "delta"]$scale
    alpha <- parmatrix[group == "alpha"]$value * parmatrix[group == "alpha"]$scale
    gamma <- parmatrix[group == "gamma"]$value * parmatrix[group == "gamma"]$scale
    eta <- parmatrix[group == "eta"]$value * parmatrix[group == "eta"]$scale
    kappa <- fgarch_moment_v(distribution = env$distribution, gamma = gamma, eta = eta, delta = delta,
                             skew = skew, shape = shape, lambda = lambda)
    alpha <- sum(alpha * kappa)
    p <- beta + alpha
    return(p)
}


.persistence_gjrgarch <- function(pars, env)
{
    parameter <- group <- value <- NULL
    parmatrix <- env$parmatrix
    parmatrix[estimate == 1, value := pars]
    alpha <- sum(parmatrix[group == "alpha"]$value * parmatrix[group == "alpha"]$scale)
    beta <- sum(parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale)
    skew <- parmatrix[parameter == "skew"]$value * parmatrix[parameter == "skew"]$scale
    shape <- parmatrix[parameter == "shape"]$value * parmatrix[parameter == "shape"]$scale
    lambda <- parmatrix[parameter == "lambda"]$value * parmatrix[parameter == "lambda"]$scale
    kappa <- gjrgarch_moment(env$distribution, skew, shape, lambda)
    gamma <- parmatrix[group == "gamma"]$value * parmatrix[group == "gamma"]$scale
    p <- alpha + beta + sum(gamma * kappa)
    return(p)
}


.persistence_cgarch <- function(pars, env)
{
    parameter <- group <- value <- NULL
    parmatrix <- env$parmatrix
    parmatrix[estimate == 1, value := pars]
    alpha <- parmatrix[group == "alpha"]$value * parmatrix[group == "alpha"]$scale
    beta <- parmatrix[group == "beta"]$value * parmatrix[group == "beta"]$scale
    rho <- parmatrix[group == "rho"]$value * parmatrix[group == "rho"]$scale
    transitory <- sum(alpha) + sum(beta)
    # [transitory, permanent]
    permanent <- rho
    return(c("permanent" = permanent, "transitory" = transitory))
}

