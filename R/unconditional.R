.unconditional_garch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$constant_variance * (1 - p)
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    return(unconditional_variance)
}

.unconditional_egarch <- function(object) {
    if (max(object$spec$model$order) > 1) {
        return(.unconditional_egarch_simulation(object))
    } else {
        return(.unconditional_egarch_analytical(object))
    }
}

.unconditional_egarch_analytical <- function(object) {
    parameter <- group <- NULL
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[parameter == "alpha1"]$value
    beta <- object$parmatrix[parameter == "beta1"]$value
    gamma <- object$parmatrix[parameter == "gamma1"]$value
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    kappa <- egarch_moment(distribution = object$spec$distribution, skew = skew, shape = shape, lambda = lambda)
    numerator <- omega
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    term1 <- exp(numerator/(1 - beta))
    prod_approx <- .egarch_moment_prod_approx(alpha = alpha, gamma = gamma,
                                              beta = beta, kappa = kappa,
                                              skew = skew, shape = shape,
                                              lambda = lambda,
                                              distribution = object$spec$distribution, n = 1000)

    moment <- term1 * prod(prod_approx)
    return(moment)

}

.unconditional_egarch_simulation <- function(object) {
    spec <- object$spec
    spec$parmatrix <- copy(object$parmatrix)
    sim <- simulate(spec, nsim = 200, h = 10000, burn = 2000)
    return(mean(rowMeans(sim$sigma^2, na.rm = TRUE), trim = 0.01))
}

.unconditional_aparch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    delta <- object$parmatrix[parameter == "delta"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    unconditional_variance <- unconditional_variance^(2/delta)
    return(unconditional_variance)
}

.unconditional_gjrgarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    return(unconditional_variance)
}

.unconditional_fgarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    delta <- object$parmatrix[parameter == "delta"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    unconditional_variance <- unconditional_variance^(2/delta)
    return(unconditional_variance)
}

.unconditional_cgarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    # the permanent component persistence
    p <- persistence(object)[1]
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    return(unconditional_variance)
}



.egarch_moment_prod_approx <- function(alpha, gamma, beta, kappa, skew, shape, lambda, distribution = "norm", n = 1000, bounds = c(-Inf, Inf))
{
    v <- sapply(1:n, function(i) {
        f <- function(x){
            terma <- exp((beta^(i - 1)) * (alpha * x + gamma * (abs(x) - kappa)))
            terma[terma > 1e4] <- 1e4
            terma[!is.finite(terma)] <- 1e4
            terma[terma < 1e-12] <- 0
            termb <- ddist(distribution, x = x, mu = 0, sigma = 1, skew = skew, shape = shape, lambda = lambda)
            termb[termb < 1e-6] <- 0
            #termb[termb > 1e4] <- 1e4
            return(terma * termb)
        }
        integrate(f, bounds[1], bounds[2], stop.on.error = FALSE)$value
    })
    return(v)
}

